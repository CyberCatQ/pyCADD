import os
from typing import List, Literal, Tuple, Union

import numpy as np

from pyCADD.utils.common import File
from pyCADD.utils.tool import NUM_PARALLEL, multiprocessing_run

from . import Structure, StructureReader, logger
from .common import DockResultFile, GridFile, MaestroFile
from .config import DataConfig
from .core import dock, grid_generate, keep_single_chain, minimize
from .data import extract_docking_data


def _write_struc(
    structure_index: int,
    structure: Structure,
    /,
    file_prefix: str = "structure",
    overwrite: bool = False,
    save_dir: str = None,
) -> MaestroFile:
    """Write a structure to a file with specified prefix and index.

    Args:
        write_info (list): a tuple containing the index and the structure. e.g. (0, <Structure>)
        file_prefix (str, optional): file prefix for generated structure file. Defaults to 'structure'.
        overwrite (bool, optional): Whether to overwrite existed file. If False, the file will be skipped. Defaults to False.
        save_dir (str, optional): directory to save the structure file. Defaults to None.

    Returns:
        str: generated structure file
    """
    save_dir = save_dir if save_dir is not None else os.getcwd()
    st_name = f"{file_prefix}-{structure_index}"
    structure.property["i_user_StructureIndex"] = structure_index
    structure.property["s_user_StructureName"] = st_name
    st_path = os.path.join(save_dir, f"{st_name}.maegz")
    if os.path.exists(st_path) and not overwrite:
        logger.debug(f"File exists: {st_path}")
        return st_path
    structure.write(st_path)
    return MaestroFile(st_path)


def split_structure(
    multi_structure_file: Union[str, File],
    save_dir: str = None,
    overwrite: bool = False,
    cpu_num: int = None,
) -> List[MaestroFile]:
    """Split a multi-structure file to single structure files.

    Args:
        multi_structure_file (Union[str, File]): path or File object of structure file containing multiple structures.
        save_dir (str, optional): directory to save the splited structure files. Defaults to None.
        overwrite (bool, optional): Whether to overwrite the existed file. Defaults to False.
        cpu_num (int, optional): cpu core number used to split structures. Defaults to 3/4 of available cores.

    Returns:
        List[MaestroFile]: list of splited structure file.
    """
    save_dir = os.path.abspath(save_dir) if save_dir is not None else multi_structure_file.file_dir
    if isinstance(multi_structure_file, str):
        multi_structure_file = File(multi_structure_file)
    os.makedirs(save_dir, exist_ok=True)
    logger.debug(f"Prepare to split structure file {multi_structure_file.file_name} to {save_dir}")

    if isinstance(multi_structure_file, str):
        multi_structure_file = File(multi_structure_file)

    structures = [st for st in StructureReader(multi_structure_file.file_path)]
    logger.debug(f"Found {len(structures)} structures in {multi_structure_file.file_name}")

    cpu_num = cpu_num if cpu_num is not None else NUM_PARALLEL

    return multiprocessing_run(
        _write_struc,
        zip(range(1, len(structures) + 1), structures),
        job_name="Split Structures",
        num_parallel=cpu_num,
        total_task_num=len(structures),
        file_prefix=multi_structure_file.file_prefix,
        save_dir=save_dir,
        overwrite=overwrite,
    )


def multi_keep_single_chain(
    structure_files: List[Union[str, MaestroFile]],
    ligand_name_list: List[str] = None,
    chain_id_list: List[str] = None,
    save_dir: str = None,
    overwrite: bool = False,
    cpu_num: int = None,
) -> List[MaestroFile]:
    """Keep the single chain for multiple structures.

    Args:
        structure_files (List[Union[str, MaestroFile]]): Structure files to be processed.
        ligand_name_list (List[str], optional): Ligand name list for each structure. Only the chain where the ligand is located will be retained.\
            Should have the same length as structure_files. Defaults to None.
        chain_id_list (List[str], optional): Chain ID list for each structure to keep. Should have the same length as structure_files. \
            Will be ignored if ligand_name_list is provided. Defaults to None.
        save_dir (str, optional): directory to save the structure files. Defaults to None.
        overwrite (bool, optional): Whether to overwrite the existed file. Defaults to False.
        cpu_num (int, optional): cpu core number used to split structures. Defaults to 3/4 of available cores.

    Raises:
        ValueError: Single-chain structure cannot be kept without ligand name.
        ValueError: The number of ligand name list is not equal to the number of structure files.

    Returns:
        List[MaestroFile]: list of structure files with single chain.
    """
    save_dir = os.path.abspath(save_dir) if save_dir is not None else os.getcwd()
    os.makedirs(save_dir, exist_ok=True)
    cpu_num = cpu_num if cpu_num is not None else NUM_PARALLEL
    logger.debug(
        f"Prepare to split and keep the single chain for {len(structure_files)} structures to {save_dir}"
    )

    structure_files = [MaestroFile(st) for st in structure_files]
    if not ligand_name_list and not chain_id_list:
        ligand_name_list = []
        for structure in structure_files:
            ligand_name = structure.metadata.ligand_name
            if not ligand_name:
                raise ValueError(
                    f"Single-chain structure cannot be kept without ligand name: {structure.file_path}"
                )
            ligand_name_list.append(ligand_name)
        if len(ligand_name_list) != len(structure_files):
            raise ValueError(
                f"The number of ligand name list ({len(ligand_name_list)}) is not equal \
                    to the number of structure files ({len(structure_files)})."
            )
        iterable = zip(structure_files, ligand_name_list)
    elif chain_id_list:
        iterable = zip(structure_files, [None] * len(structure_files), chain_id_list)
    else:
        iterable = zip(structure_files, ligand_name_list)

    return multiprocessing_run(
        keep_single_chain,
        job_name="Keep Single Chain",
        iterable=iterable,
        num_parallel=cpu_num,
        total_task_num=len(structure_files),
        save_dir=save_dir,
        overwrite=overwrite,
    )


def multi_minimize(
    structure_files: List[Union[str, File]],
    ph: float = 7.4,
    force_field: Literal["OPLS4", "OPLS3e", "OPLS3", "OPLS_2005"] = "OPLS4",
    fill_side_chain: bool = True,
    add_missing_loop: bool = True,
    del_water: bool = True,
    watdist: float = 5.0,
    rmsd_cutoff: float = 0.3,
    save_dir: str = None,
    overwrite: bool = False,
    cpu_num: int = None,
) -> List[MaestroFile]:
    """Minimize multiple structures.

    Args:
        structure_files (List[Union[str, File]]): list of structure file paths or File objects.
        ph (float, optional): pH value to calculate protonation states. Defaults to 7.4.
        force_field (str): force field to use. Defaults to 'OPLS4'.
        fill_side_chain (bool, optional): whether to fill side chain. Defaults to True.
        add_missing_loop (bool, optional): whether to add missing loop. Defaults to True.
        del_water (bool, optional): whether to delete water molecules. Defaults to True.
        watdist (float, optional): how far from the ligand to delete water molecules. Set to 0.0 to delete all water molecules. Defaults to 5.0.
        rmsd_cutoff (float, optional): RMSD cutoff for minimization. Defaults to 0.3.
        save_dir (str, optional): directory to save the minimized structures. Defaults to None.
        overwrite (bool, optional): Whether to overwrite the existed file. Defaults to False.
        cpu_num (int, optional): cpu core number used to minimize structures. Defaults to 3/4 of available cores.

    Returns:
        List[MaestroFile]: list of minimized structure file.
    """
    save_dir = os.path.abspath(save_dir) if save_dir is not None else os.getcwd()
    os.makedirs(save_dir, exist_ok=True)
    logger.debug(
        f"Prepare to optimize and minimize {len(structure_files)} structures at {save_dir}"
    )

    if isinstance(structure_files[0], str):
        structure_files = [File(st) for st in structure_files]

    cpu_num = cpu_num if cpu_num is not None else NUM_PARALLEL

    return multiprocessing_run(
        minimize,
        structure_files,
        job_name="Minimize",
        num_parallel=cpu_num,
        total_task_num=len(structure_files),
        ph=ph,
        force_field=force_field,
        fill_side_chain=fill_side_chain,
        add_missing_loop=add_missing_loop,
        del_water=del_water,
        watdist=watdist,
        rmsd_cutoff=rmsd_cutoff,
        save_dir=save_dir,
        overwrite=overwrite,
    )


def multi_grid_generate(
    structure_files: List[Union[str, File]],
    box_center_list: List[Tuple[float, float, float]] = None,
    box_center_molnum_list: List[int] = None,
    box_size_list: List[float] = None,
    force_field: Literal["OPLS4", "OPLS3e", "OPLS3", "OPLS_2005"] = "OPLS4",
    save_dir: str = None,
    overwrite: bool = False,
    cpu_num: int = None,
) -> List[GridFile]:
    """Generate multiple grid files.

    Args:
        structure_files (List[Union[str, File]]): list of structure file paths or File objects.
        box_center_list (List[Tuple[float, float, float]], optional): list of box centers. Should be the same length as structure_files. Defaults to None.
        box_center_molnum_list (List[int], optional): list of molecule numbers for box centers. Should be the same length as structure_files. Defaults to None.
        box_size_list: list of box sizes. Should be the same length as structure_files. Defaults to [20, 20, ...] with the length of structure_files.
        force_field (str, optional): force field to use. Defaults to 'OPLS4'.
        save_dir (str, optional): directory to save the grid files. Defaults to None.
        overwrite (bool, optional): Whether to overwrite the existed file. Defaults to False.
        cpu_num (int, optional): cpu core number used to generate grid files. Defaults to 3/4 of available cores.

    Returns:
        List[GridFile]: list of grid file.
    """
    cpu_num = cpu_num if cpu_num is not None else NUM_PARALLEL
    save_dir = os.path.abspath(save_dir) if save_dir is not None else os.getcwd()
    structure_length = len(structure_files)

    if box_size_list is None:
        box_size_list = [20] * structure_length
    elif len(box_size_list) != structure_length:
        raise ValueError(f"box_size_list should have the same length as structure_files.")

    zip_args = None
    if box_center_list is not None:
        if np.array(box_center_list).shape != (structure_length, 3):
            raise ValueError(
                f"box_center_list should have the same length as structure_files and each element should be a tuple of 3 floats."
            )
        zip_args = zip(
            structure_files, box_center_list, [None for i in range(structure_length)], box_size_list
        )
    elif box_center_molnum_list is not None:
        if len(box_center_molnum_list) != structure_length:
            raise ValueError(
                f"box_center_molnum_list should have the same length as structure_files."
            )
        zip_args = zip(
            structure_files,
            [None for i in range(structure_length)],
            box_center_molnum_list,
            box_size_list,
        )
    else:
        raise ValueError("Either box_center_list or box_center_molnum_list should be provided.")

    if isinstance(structure_files[0], str):
        structure_files = [File(st) for st in structure_files]

    os.makedirs(save_dir, exist_ok=True)
    logger.debug(
        f"Prepare to generate grid files for {len(structure_files)} structures at {save_dir}"
    )

    return multiprocessing_run(
        grid_generate,
        zip_args,
        job_name="Generate Grids",
        num_parallel=cpu_num,
        total_task_num=len(structure_files),
        force_field=force_field,
        save_dir=save_dir,
        overwrite=overwrite,
    )


def get_docking_pairs(
    grid_files: List[Union[str, File]], ligand_files: List[Union[str, File]]
) -> List[Tuple[Union[str, File], Union[str, File]]]:
    """Get all possible docking pairs from grid and ligand files, \
        which will establish a mapping relationship between each provided receptor and each provided ligand.

    Args:
        grid_files (List[Union[str, File]]): Grid files.
        ligand_files (List[Union[str, File]]): Ligand files.

    Returns:
        tuple[Union[str, File], Union[str, File]: a tuple containing grid and ligand file. 
    """
    return [(grid, ligand) for grid in grid_files for ligand in ligand_files]


def multi_dock(
    docking_pairs: List[Tuple[GridFile, File]],
    force_field: Literal["OPLS4", "OPLS3e", "OPLS3", "OPLS_2005"] = "OPLS4",
    precision: Literal["SP", "XP", "HTVS"] = "SP",
    calc_rmsd: bool = False,
    include_receptor: bool = False,
    save_dir: str = None,
    overwrite: bool = False,
    cpu_num: int = None,
) -> List[DockResultFile]:
    """Dock multiple ligands to multiple grids.

    Args:
        docking_mapping (List[Tuple[GridFile, File]]): list of tuples containing grid and ligand file paths or File objects.\
            e.g. [(grid1, ligand1), (grid2, ligand2)]
        force_field (str, optional): force field to use. Defaults to 'OPLS4'.
        precision (str, optional): docking precision. Defaults to 'SP'.
        calc_rmsd (bool, optional): Whether to calculate RMSD with co-crystal ligand. \
            If True, grid file must be generated from a complex. Defaults to False.
        include_receptor (bool, optional): Whether to include receptor structure in the output file. Defaults to False.
        save_dir (str, optional): directory to save the docked results. \
            Results will be further saved in separate directories named PDBIDs. Defaults to None.
        overwrite (bool, optional): Whether to overwrite the existed file. Defaults to False.
        cpu_num (int, optional): cpu core number used to dock ligands. Defaults to 3/4 of available cores.

    Returns:
        List[DockResultFile]: list of docked result file.
    """
    cpu_num = cpu_num if cpu_num is not None else NUM_PARALLEL
    save_dir = os.path.abspath(save_dir) if save_dir is not None else os.getcwd()
    os.makedirs(save_dir, exist_ok=True)
    logger.debug(f"Prepare to dock {len(docking_pairs)} pairs at {save_dir}")

    return multiprocessing_run(
        dock,
        docking_pairs,
        job_name="Dock",
        num_parallel=cpu_num,
        force_field=force_field,
        precision=precision,
        calc_rmsd=calc_rmsd,
        include_receptor=include_receptor,
        save_dir=save_dir,
        overwrite=overwrite,
    )


def multi_extract_data(
    dock_result_files: List[Union[str, File]], data_config: DataConfig = None, cpu_num: int = None
) -> List[dict]:
    """Extract data from multiple docked result files.

    Args:
        dock_result_files (List[Union[str, File]]): list of docked result file paths or File objects.
        extract_type (str, optional): data type to extract. Defaults to 'all'.
        save_dir (str, optional): directory to save the extracted data. Defaults to None.
        cpu_num (int, optional): cpu core number used to extract data. Defaults to 3/4 of available cores.

    Returns:
        List[dict]: list of extracted data dict, which can be converted to DataFrame directly.
    """
    cpu_num = cpu_num if cpu_num is not None else NUM_PARALLEL
    logger.debug(f"Prepare to extract data from {len(dock_result_files)} docked result files")

    nested_result = multiprocessing_run(
        extract_docking_data,
        dock_result_files,
        job_name="Extract Data",
        num_parallel=cpu_num,
        total_task_num=len(dock_result_files),
        data_config=data_config,
    )
    return [data for data_list in nested_result for data in data_list]
