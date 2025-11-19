import os
import shutil
from typing import Literal, Tuple, Union

from pyCADD.Dock.schrodinger.common import DockResultFile, GridFile, MaestroFile
from pyCADD.Dock.schrodinger.config import FORCE_FILED_DICT
from pyCADD.Dock.schrodinger.utils import convert_format, launch
from pyCADD.utils.common import ChDir
from pyCADD.utils.tool import shell_run

from . import Structure, evaluate_asl, logger


def split_complex(
    file: Union[MaestroFile, str],
    ligand_atom_indexes: list = None,
    ligand_asl: str = None,
    structure_index: int = 0,
) -> Tuple[Structure, Structure]:
    """Split the complex structure into receptor and ligand.

    Args:
        file (MaestroFile | str): Complex structure file object or path.
        ligand_atom_indexes (list, optional): List of ligand atom indexes. Defaults to None.
        ligand_asl (str, optional): ASL to identify ligand. Defaults to None.
        structure_index (int, optional): Index of the structure to split. Defaults to 0.

    Raises:
        ValueError: Ligand_atom_indexes and ligand_asl cannot be both None.
        ValueError: Atom indexes or ASL identified multiple molecules or no atom.

    Returns:
        Tuple[Structure, Structure]: Receptor and Ligand structures.
    """
    if isinstance(file, str):
        file = MaestroFile(path=file)
    st = file.structures[structure_index]

    atom_indexes = None
    if ligand_atom_indexes is not None:
        atom_indexes = ligand_atom_indexes
    elif ligand_asl is not None:
        atom_indexes = evaluate_asl(st, ligand_asl)
    else:
        raise ValueError("Please provide ligand_atom_indexes or ligand_asl.")

    mol_set = set(st.atom[index].molecule_number for index in atom_indexes)

    if len(mol_set) != 1:
        msg = "Identified "
        if len(mol_set) > 1:
            msg += "multiple molecules"
        elif len(mol_set) == 0:
            msg += "no atom"
        msg += " with the ligand_atom_indexes or ligand_asl."
        raise ValueError(msg)

    ligand_st = st.extract(atom_indexes)
    receptor_st = st.copy()
    receptor_st.deleteAtoms(atom_indexes)
    return receptor_st, ligand_st


def minimize(
    file: Union[MaestroFile, str],
    ph: float = 7.4,
    force_field: Literal["OPLS4", "OPLS3e", "OPLS3", "OPLS_2005"] = "OPLS4",
    fill_side_chain: bool = True,
    add_missing_loop: bool = True,
    del_water: bool = True,
    watdist: float = 5.0,
    rmsd_cutoff: float = 0.3,
    save_dir: str = None,
    overwrite: bool = False,
) -> MaestroFile:
    """Prepare the structure and run minimization using Schrodinger's prepwizard.

    Args:
        file (Union[MaestroFile, str]): file path or MaestroFile object to be prepared and minimized.
        ph (float, optional): pH value to calculate protonation states. Defaults to 7.4.
        force_field (str): force field to use. Defaults to 'OPLS4'.
        fill_side_chain (bool, optional): whether to fill side chain. Defaults to True.
        add_missing_loop (bool, optional): whether to add missing loop. Defaults to True.
        del_water (bool, optional): whether to delete water molecules. Defaults to True.
        watdist (float, optional): how far from the ligand to delete water molecules. Set to 0.0 to delete all water molecules. Defaults to 5.0.
        rmsd_cutoff (float, optional): RMSD cutoff for minimization. Defaults to 0.3.
        save_dir (str, optional): directory to save the results. Defaults to None.
        overwrite (bool, optional): whether to overwrite existing result files. Defaults to False.

    Raises:
        RuntimeError: raise if the minimization process failed.

    Returns:
        MaestroFile: minimized structure file
    """
    if isinstance(file, str):
        file = MaestroFile(path=file)
    save_dir = os.path.abspath(save_dir) if save_dir else os.getcwd()
    os.makedirs(save_dir, exist_ok=True)
    metadata = file.metadata.copy()
    metadata.action = "minimized"
    file_prefix = metadata.generate_file_name(["pdbid", "ligand_name", "action"])

    force_field = FORCE_FILED_DICT[force_field]
    logger.debug(f"Prepare to prepare and minimize file: {file.file_path}")
    minimized_file = os.path.join(save_dir, f"{file_prefix}.mae")
    if not overwrite and os.path.exists(minimized_file):
        logger.debug(f"File already existed: {minimized_file}")
        return MaestroFile(minimized_file, metadata=metadata)

    with ChDir(save_dir):
        job_name = f"{file.file_prefix}-Minimize"
        prepwizard_command = f"prepwizard -f {force_field} -r {rmsd_cutoff} -propka_pH {ph} -disulfides -s -JOBNAME {job_name}"
        if fill_side_chain:
            prepwizard_command += " -fillsidechains"
        if add_missing_loop:
            prepwizard_command += " -fillloops"
        if del_water:
            prepwizard_command += f" -watdist {watdist}"
        else:
            prepwizard_command += " -keepfarwat"
        _minimized_file = f"_{file_prefix}.mae"
        prepwizard_command += f" {file.file_path} {_minimized_file}"

        launch(prepwizard_command)

        try:
            shutil.move(_minimized_file, minimized_file)
        except FileNotFoundError:
            raise RuntimeError(f"Minimization Process Failed: {file.file_prefix}")
        else:
            logger.debug(f"Minimized file saved: {minimized_file}")
            return MaestroFile(minimized_file, metadata=metadata)


def grid_generate(
    file: Union[MaestroFile, str],
    box_center: Tuple[float, float, float] = None,
    box_center_molnum: int = None,
    box_size: int = 20,
    force_field: Literal["OPLS4", "OPLS3e", "OPLS3", "OPLS_2005"] = "OPLS4",
    save_dir: str = None,
    overwrite: bool = False,
) -> GridFile:
    """Generate grid file for Glide docking.

    Args:
        file (Union[MaestroFile, str]): structure to generate grid file.
        box_center (tuple[float, float, float], optional): center XYZ of the grid box. Defaults to None.
        box_center_molnum (int, optional): the molecule number of the molecule that is set as the center.\
            This molecule will be removed during the grid box generation process,\
            and its centroid will be used as the center of the box.\
            Will be ignored when box_center is set. Defaults to None.
        box_size (int, optional): box size of grid. Defaults to 20.
        force_field (str, optional): force field to use. Defaults to 'OPLS4'.
        save_dir (str, optional): directory to save the results. Defaults to None.
        overwrite (bool, optional): whether to overwrite existing files. Defaults to False.

    Raises:
        RuntimeError: raise if the grid generation process failed.

    Returns:
        GridFile: grid file
    """
    save_dir = os.path.abspath(save_dir) if save_dir is not None else os.getcwd()
    os.makedirs(save_dir, exist_ok=True)
    if isinstance(file, str):
        file = MaestroFile(path=file)
    if file.file_ext not in ["mae", "maegz"]:
        file = MaestroFile(
            convert_format(file.file_path, "mae", save_dir=file.file_dir), metadata=file.metadata
        )
    metadata = file.metadata.copy()
    metadata.action = "glide-grid"
    file_prefix = metadata.generate_file_name(["pdbid", "ligand_name", "action"])
    grid_file = os.path.join(save_dir, f"{file_prefix}.zip")
    logger.debug(f"Prepare to generate grid file: {grid_file}")

    if os.path.exists(grid_file) and not overwrite:
        logger.debug(f"File already existed: {grid_file}")
        return GridFile(grid_file, metadata=metadata)

    with ChDir(save_dir):
        input_file = f"{file_prefix}.in"
        job_name = file_prefix.replace("_", "-")
        outsize = box_size + 10

        glide_grid_config = [
            f"FORCEFIELD {force_field}\n",
            f"GRIDFILE {grid_file}\n",
            "INNERBOX 10,10,10\n",
            f"OUTERBOX {outsize},{outsize},{outsize} \n",
            f"RECEP_FILE {file.file_path}\n",
        ]
        if box_center is not None:
            try:
                coord_x, coord_y, coord_z = box_center
            except ValueError:
                raise ValueError("box_center should be a tuple of 3 floats.")
            glide_grid_config += [f"GRID_CENTER {coord_x},{coord_y},{coord_z}\n"]
            logger.debug(f"Grid box center is set to {box_center}")
        elif box_center_molnum is not None:
            glide_grid_config += [f"LIGAND_MOLECULE {box_center_molnum}\n"]
            logger.debug(f"Grid box center is set to Molecule {box_center_molnum}")
        else:
            raise ValueError("Either box_center or box_center_molnum should be provided.")

        with open(input_file, "w") as f:
            f.writelines(glide_grid_config)

        shell_run(f"glide {input_file} -JOBNAME {job_name} -WAIT -NOJOBID > {job_name}.log 2>&1")
        if not os.path.exists(grid_file):
            raise RuntimeError(f"{grid_file} Generation Failed.")
        logger.debug(f"Grid file generated: {grid_file}")
    return GridFile(grid_file, metadata=metadata)


def dock(
    grid_file: Union[GridFile, str],
    ligand_file: Union[MaestroFile, str],
    force_field: Literal["OPLS4", "OPLS3e", "OPLS3", "OPLS_2005"] = "OPLS4",
    precision: Literal["SP", "XP", "HTVS"] = "SP",
    calc_rmsd: bool = False,
    include_receptor: bool = False,
    save_dir: str = None,
    overwrite: bool = False,
) -> DockResultFile:
    """Perform Glide ligand docking

    Args:
        grid_file (Union[GridFile, str]): grid file object or path.
        ligand_file (Union[MaestroFile, str]): Ligand file object or path.
        force_field (str, optional): force field to use. Defaults to 'OPLS4'.
        precision (str, optional): docking precision. Defaults to 'SP'.
        calc_rmsd (bool, optional): Whether to calculate RMSD with co-crystal ligand. If True, grid file must be generated from a complex. Defaults to False.
        include_receptor (bool, optional): Whether to include receptor structure in the output file. Defaults to False.
        save_dir (str, optional): directory to save the results. Results will be further saved in directory named PDBID. \
            Using save_dir directly if PDBID is not provided. Defaults to current working directory.
        overwrite (bool, optional): whether to overwrite existing files. Defaults to False.

    Raises:
        RuntimeError: raise if the docking process failed.

    Returns:
        DockResultFile: docking result file.
    """

    if isinstance(grid_file, str):
        grid_file = GridFile(grid_file)
    if isinstance(ligand_file, str):
        ligand_file = MaestroFile(ligand_file)

    save_dir = os.path.abspath(save_dir) if save_dir is not None else os.getcwd()
    save_dir = (
        os.path.join(save_dir, grid_file.metadata.pdbid) if grid_file.metadata.pdbid else save_dir
    )
    os.makedirs(save_dir, exist_ok=True)

    if ligand_file.file_ext not in ["mae", "maegz"]:
        coverted_ligand_file = convert_format(
            ligand_file.file_path,
            f"{ligand_file.file_prefix}.maegz",
            save_dir=save_dir,
            overwrite=True,
        )
        ligand_file = MaestroFile(coverted_ligand_file, metadata=ligand_file.metadata)

    metadata = grid_file.metadata.copy()
    metadata.docking_ligand_name = ligand_file.file_prefix.replace("_", "-")
    metadata.precision = precision
    metadata.set("calc_rmsd", bool(calc_rmsd))
    metadata.set("include_receptor", bool(include_receptor))
    metadata.action = "glide-dock"

    file_prefix = metadata.generate_file_name(
        ["pdbid", "ligand_name", "action", "docking_ligand_name", "precision"]
    )
    dock_result_file = os.path.join(save_dir, f"{file_prefix}.maegz")

    if os.path.exists(dock_result_file) and not overwrite:
        logger.debug(f"File already existed: {dock_result_file}")
        return DockResultFile(dock_result_file, metadata=metadata)

    with ChDir(save_dir):
        input_file = f"{file_prefix}.in"
        job_name = file_prefix.replace("_", "-")
        output_file = f"{job_name}_pv.maegz" if include_receptor else f"{job_name}_lib.maegz"

        glide_dock_config = [
            f"GRIDFILE {grid_file.file_path}\n",
            f"LIGANDFILE {ligand_file.file_path}\n",
            f"PRECISION {precision}\n",
            f"FORCEFIELD {force_field}\n",
        ]
        if not include_receptor:
            glide_dock_config.append("POSE_OUTTYPE ligandlib\n")
        if calc_rmsd:
            glide_dock_config.append("CALC_INPUT_RMS True\n")
        if precision == "XP":
            glide_dock_config.append("WRITE_XP_DESC False\n")
            glide_dock_config.append("POSTDOCK_XP_DELE 0.5\n")

        with open(input_file, "w") as f:
            f.writelines(glide_dock_config)

        shell_run(f"glide {input_file} -JOBNAME {job_name} -WAIT -NOJOBID > {job_name}.log 2>&1")

        try:
            shutil.move(output_file, dock_result_file)
        except FileNotFoundError:
            logger.debug(f"Docking Failed: {file_prefix}")
            return DockResultFile(dock_result_file, metadata=metadata, exist=False)

        logger.debug(f"Docking result file saved: {dock_result_file}")

    return DockResultFile(dock_result_file, metadata=metadata)


def keep_single_chain(
    structure_file: Union[MaestroFile, str],
    by_resname: str = None,
    chain_id: str = None,
    save_dir: str = None,
    overwrite: bool = False,
) -> MaestroFile:
    """Keep single chain from a structure file.

    Args:
        structure_file (Union[MaestroFile, str]): Structure file object or path.
        by_resname (str, optional): Residue or ligand name from the chain to be kept. Defaults to None.
        chain_id (str, optional): Chain ID to be kept. Defaults to None.
        save_dir (str, optional): Directory to save the results. Defaults to None.
        overwrite (bool, optional): Whether to overwrite existing files. Defaults to False.

    Raises:
        ValueError: Either chain_id or by_resname should be specified.

    Returns:
        MaestroFile: The structure file with single chain.
    """
    save_dir = save_dir if save_dir is not None else os.getcwd()
    if isinstance(structure_file, str):
        structure_file = MaestroFile(structure_file)
    metadata = structure_file.metadata.copy()
    metadata.action = "keep-single-chain"

    if chain_id is None:
        if by_resname is not None:
            ligand = structure_file.find_ligands(specified_names=[by_resname])
            if len(ligand) > 1:
                msg = f"{len(ligand)} ligands named {by_resname} found in the structure file, and the first one will be kept.\n"
                msg += f"Keep the chain structure of ligand {ligand[0]}"
                logger.warning(msg)
            chain_id = ligand[0].chain
            metadata.ligand_name = by_resname
        else:
            raise ValueError("Please specify either chain_id or by_resname")
    metadata.pdbid = f"{structure_file.file_prefix}-chain-{chain_id}"
    file_prefix = metadata.generate_file_name(["pdbid", "ligand_name"])
    file_suffix = structure_file.file_suffix
    chain_structure_file = os.path.join(save_dir, f"{file_prefix}.{file_suffix}")
    if not overwrite and os.path.exists(chain_structure_file):
        logger.debug(f"File already existed: {chain_structure_file}")
        return MaestroFile(chain_structure_file, metadata=metadata)
    chain_structure = structure_file.get_chain_structure(chain_id=chain_id)
    chain_structure.write(chain_structure_file)
    return MaestroFile(chain_structure_file, metadata=metadata)
