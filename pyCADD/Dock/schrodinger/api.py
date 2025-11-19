import os
from datetime import datetime
from typing import List, Literal, Tuple, Union

import pandas as pd

from pyCADD.Dock.common import EnsembleInputFile
from pyCADD.Dock.schrodinger.common import DockResultFile, GridFile, MaestroFile, MetaData
from pyCADD.Dock.schrodinger.core import dock, grid_generate, minimize, split_complex
from pyCADD.Dock.schrodinger.ensemble import (
    get_docking_pairs,
    multi_dock,
    multi_extract_data,
    multi_grid_generate,
    multi_keep_single_chain,
    multi_minimize,
    split_structure,
)
from pyCADD.utils.common import ChDir
from pyCADD.utils.tool import download_pdb, download_pdb_list
from . import logger


DATE = datetime.now().strftime("%Y%m%d%H%M")


class DockControl:
    """Ligand Dock Control Class"""

    def __init__(self, pdbid: str = None, save_path: str = None) -> None:
        """Ligand Dock Control Class

        Args:
            pdbid (str, optional): PDB ID from RCSB. Defaults to None.
            save_path (str, optional): directory to save the result files. Defaults to None.

        Raises:
            ValueError: either protein_file or pdbid must be provided
        """
        self.save_path = os.path.abspath(save_path) if save_path is not None else os.getcwd()
        os.makedirs(self.save_path, exist_ok=True)

        self.protein_file = self._download_pdb_file(pdbid) if pdbid else None
        self.minimized_file = None
        self.grid_file = None
        self.dock_result_file = None

    def _download_pdb_file(self, pdbid: str) -> str:
        """Download PDB file from RCSB PDB.

        Args:
            pdbid (str): PDB ID

        Raises:
            RuntimeError: Failed to download pdb file from RCSB PDB.

        Returns:
            str: pdb file path
        """
        with ChDir(self.save_path):
            download_pdb(pdbid, save_dir="pdb")
            pdb_file = os.path.join(self.save_path, "pdb", f"{pdbid}.pdb")
            if not os.path.exists(pdb_file):
                raise RuntimeError(f"Failed to download pdb file from RCSB PDB: {pdbid}")
        return pdb_file

    def split(
        self,
        structure_file: str = None,
        ligand_resname: str = None,
        ligand_atom_indexes: list = None,
        ligand_asl: str = None,
    ) -> Tuple[str, str]:
        """Split structure to protein and ligand.

        Args:
            structure_file (str, optional): structure file. Defaults to None.
            ligand_resname (str, optional): specify the ligand name to be splited. Defaults to None.
            ligand_atom_indexes (list, optional): atom indexes to define the ligand molecule. Defaults to None.
            ligand_asl (str, optional): ASL to define the ligand molecule. Defaults to None.

        Returns:
            Tuple[str, str]: protein_file_path, ligand_file_path
        """
        if structure_file is None:
            structure_file = (
                self.protein_file if self.minimized_file is None else self.minimized_file
            )
        if isinstance(structure_file, str):
            structure_file = MaestroFile(structure_file)
        if not ligand_atom_indexes and not ligand_asl:
            lig_name = [ligand_resname] if ligand_resname is not None else None
            ligand = structure_file.find_ligands(specified_names=lig_name)[0]
            ligand_atom_indexes = ligand.atom_indexes
        with ChDir(self.save_path):
            pro_st, lig_st = split_complex(
                structure_file, ligand_atom_indexes=ligand_atom_indexes, ligand_asl=ligand_asl
            )
            os.makedirs("protein", exist_ok=True)
            os.makedirs("ligand", exist_ok=True)
            pro_file = os.path.join("protein", f"{structure_file.file_prefix}_pro.pdb")
            lig_file = os.path.join("ligand", f"{structure_file.file_prefix}_lig.sdf")
            pro_st.write(pro_file)
            lig_st.write(lig_file)
        return pro_file, lig_file

    def minimize(
        self,
        structure_file: Union[str, MaestroFile] = None,
        ph: float = 7.4,
        force_field: Literal["OPLS4", "OPLS3e", "OPLS3", "OPLS_2005"] = "OPLS4",
        fill_side_chain: bool = True,
        add_missing_loop: bool = True,
        del_water: bool = True,
        watdist: float = 5.0,
        rmsd_cutoff: float = 0.3,
        overwrite: bool = False,
    ) -> MaestroFile:
        """Minimize the protein structure.

        Args:
            structure_file (Union[MaestroFile, str], optional): structure to minimize. If not specified, use the initial protein structure. Defaults to None.
            ph (float, optional): pH value to calculate protonation states. Defaults to 7.4.
            force_field (str): force field to use. Defaults to 'OPLS4'.
            fill_side_chain (bool, optional): whether to fill side chain. Defaults to True.
            add_missing_loop (bool, optional): whether to add missing loop. Defaults to True.
            del_water (bool, optional): whether to delete water molecules. Defaults to True.
            watdist (float, optional): how far from the ligand to delete water molecules. Set to 0.0 to delete all water molecules. Defaults to 5.0.
            rmsd_cutoff (float, optional): RMSD cutoff for minimization. Defaults to 0.3.
            overwrite (bool, optional): whether to overwrite existing result files. Defaults to False.

        Returns:
            MaestroFile: Minimized structure file.
        """
        with ChDir(self.save_path):
            self.minimized_file = minimize(
                structure_file if structure_file is not None else self.protein_file,
                ph=ph,
                force_field=force_field,
                fill_side_chain=fill_side_chain,
                add_missing_loop=add_missing_loop,
                del_water=del_water,
                watdist=watdist,
                rmsd_cutoff=rmsd_cutoff,
                save_dir="minimize",
                overwrite=overwrite,
            )
        return self.minimized_file

    def grid_generate(
        self,
        structure_file: Union[str, MaestroFile] = None,
        box_center: Tuple[float, float, float] = None,
        box_center_molnum: int = None,
        box_center_resname: str = None,
        box_size: int = 20,
        force_field: Literal["OPLS4", "OPLS3e", "OPLS3", "OPLS_2005"] = "OPLS4",
        overwrite: bool = False,
    ) -> GridFile:
        """Generate grid file for docking.

        Args:
            structure_file (Union[MaestroFile, str]): structure to generate grid file. If not specified, use the minimized structure.
            box_center (tuple[float, float, float], optional): center XYZ of the grid box. Defaults to None.
            box_center_molnum (int, optional): the molecule number of the molecule that is set as the center.\
                This molecule will be removed during the grid box generation process,\
                and its centroid will be used as the center of the box.\
                Will be ignored when box_center is set. Defaults to None.
            box_center_resname (str, optional): the residue name of the molecule that is set as the center.\
                This molecule will be removed during the grid box generation process,\
                and its centroid will be used as the center of the box.\
                Will be ignored when box_center is set. Defaults to None.
            box_size (int, optional): box size of grid. Defaults to 20.
            force_field (str, optional): force field to use. Defaults to 'OPLS4'.
            overwrite (bool, optional): whether to overwrite existing files. Defaults to False.

        Returns:
            GridFile: generated grid file.
        """
        if structure_file is not None:
            self.minimized_file = (
                MaestroFile(structure_file) if isinstance(structure_file, str) else structure_file
            )
        if not box_center and not box_center_molnum:
            box_center_resname = [box_center_resname] if box_center_resname is not None else None
            box_center_molnum = self.minimized_file.find_ligands(box_center_resname)[0].mol_num

        with ChDir(self.save_path):
            self.grid_file = grid_generate(
                self.minimized_file,
                box_center=box_center,
                box_center_molnum=box_center_molnum,
                box_size=box_size,
                force_field=force_field,
                save_dir="grid",
                overwrite=overwrite,
            )
        return self.grid_file

    def dock(
        self,
        ligand_file: Union[str, MaestroFile],
        grid_file: GridFile = None,
        force_field: Literal["OPLS4", "OPLS3e", "OPLS3", "OPLS_2005"] = "OPLS4",
        precision: Literal["SP", "XP", "HTVS"] = "SP",
        calc_rmsd: bool = False,
        include_receptor: bool = False,
        overwrite: bool = False,
    ) -> DockResultFile:
        """Perform molecule docking.

        Args:
            ligand_file (Union[MaestroFile, str]): Ligand file object or path.
            grid_file (Union[GridFile, str]): grid file object or path. If not specified, use the grid file generated before.
            force_field (str, optional): force field to use. Defaults to 'OPLS4'.
            precision (str, optional): docking precision. Defaults to 'SP'.
            calc_rmsd (bool, optional): Whether to calculate RMSD with co-crystal ligand. If True, grid file must be generated from a complex. Defaults to False.
            include_receptor (bool, optional): Whether to include receptor structure in the output file. Defaults to False.
            overwrite (bool, optional): whether to overwrite existing files. Defaults to False.

        Returns:
            DockResultFile: Docking result file.
        """
        grid_file = grid_file if grid_file is not None else self.grid_file
        with ChDir(self.save_path):
            self.dock_result_file = dock(
                grid_file=grid_file,
                ligand_file=ligand_file,
                force_field=force_field,
                precision=precision,
                calc_rmsd=calc_rmsd,
                include_receptor=include_receptor,
                save_dir="dock",
                overwrite=overwrite,
            )
        return self.dock_result_file


class DockEnsemble:
    """
    Docking Control for Ensemble Docking
    """

    def __init__(
        self, input_file: Union[str, EnsembleInputFile], save_path: str = None, cpu_num: int = None
    ) -> None:
        """Docking Control for Ensemble Docking

        Args:
            input_file (Union[str, EnsembleInputFile]): Input file object or path for ensemble docking.
            save_path (str, optional): directory to save the result files. Defaults to None.
            cpu_num (int, optional): number of CPU cores to use in parallel. Defaults to 3 / 4 cores.
        """
        self.input_file = (
            EnsembleInputFile.parse_file(input_file) if isinstance(input_file, str) else input_file
        )
        self.pairs_list = self.input_file.get_pairs_list()
        self.pdbid_list = self.input_file.get_pdbid_list()

        self.pdb_files = []
        self.minimized_files = []
        self.grid_files = []
        self.cocrystal_ligands = []
        self.library_files = []
        self.ensemble_docking_results = []

        self.save_path = os.path.abspath(save_path) if save_path is not None else os.getcwd()
        os.makedirs(self.save_path, exist_ok=True)
        self.cpu_num = cpu_num if cpu_num is not None else os.cpu_count() // 4 * 3
        self._download_pdb()

    def _download_pdb(self, overwrite: bool = False) -> None:
        """Download PDB files from RCSB PDB

        Args:
            overwrite (bool, optional): Whether to overwrite existing files. Defaults to False.
        Raises:
            RuntimeError: Failed to download pdb file from RCSB PDB.
        """
        with ChDir(self.save_path):
            os.makedirs("pdb", exist_ok=True)
            download_pdb_list(self.pdbid_list, save_dir="pdb", overwrite=overwrite)
            for pdbid in self.pdbid_list:
                if not os.path.exists(os.path.join("pdb", f"{pdbid}.pdb")):
                    raise RuntimeError(f"Failed to download pdb file from RCSB PDB: {pdbid}")
            self.pdb_files = [
                MaestroFile(
                    path=os.path.join("pdb", f"{pdbid}.pdb"),
                    metadata=MetaData(pdbid=pdbid, ligand_name=lig_name),
                )
                for pdbid, lig_name in self.pairs_list
            ]

    def load_library(
        self, ligand_file: Union[str, MaestroFile], overwrite: bool = False
    ) -> List[MaestroFile]:
        """Load the compound/ligand library and split them into multiple single structure files for ensemble docking.

        Args:
            ligand_file (Union[str, MaestroFile]): library file object or path.
            overwrite (bool, optional): Whether to overwrite existing files. Defaults to False.

        Returns:
            List[MaestroFile]: Splitted ligands files list
        """
        if isinstance(ligand_file, str):
            ligand_file = MaestroFile(ligand_file)
        with ChDir(self.save_path):
            os.makedirs("ligands", exist_ok=True)
            self.library_files = split_structure(
                ligand_file, save_dir="ligands", overwrite=overwrite, cpu_num=self.cpu_num
            )
            logger.info(
                f"{len(self.library_files)} ligands have been loaded from {ligand_file.file_path}"
            )
        return self.library_files

    def keep_single_chain(self, overwrite: bool = False):
        """Keep the single chain in the structure file for ensemble docking.

        Args:
            overwrite (bool, optional): Whether to overwrite existing files. Defaults to False.
        """
        with ChDir(self.save_path):
            self.pdb_files = multi_keep_single_chain(
                self.pdb_files, save_dir="pdb", overwrite=overwrite, cpu_num=self.cpu_num
            )

    def minimize(
        self,
        ph: float = 7.4,
        force_field: Literal["OPLS4", "OPLS3e", "OPLS3", "OPLS_2005"] = "OPLS4",
        fill_side_chain: bool = True,
        add_missing_loop: bool = True,
        del_water: bool = True,
        watdist: float = 5.0,
        rmsd_cutoff: float = 0.3,
        overwrite: bool = False,
    ) -> List[MaestroFile]:
        """Minimize the structure file for ensemble docking.

        Args:
            ph (float, optional): pH value for the structure file. Defaults to 7.4.
            force_field (Literal["OPLS4", "OPLS3e", "OPLS3", "OPLS_2005"], optional): Force field for the structure file. Defaults to "OPLS4".
            fill_side_chain (bool, optional): Whether to fill side chain. Defaults to True.
            add_missing_loop (bool, optional): Whether to add missing loop. Defaults to True.
            del_water (bool, optional): Whether to delete water. Defaults to True.
            watdist (float, optional): Water distance cutoff. Defaults to 5.0.
            rmsd_cutoff (float, optional): RMSD cutoff. Defaults to 0.3.
            overwrite (bool, optional): Whether to overwrite existing files. Default to False.

        Returns:
            List[MaestroFile]: Minimized structure files list
        """

        with ChDir(self.save_path):
            self.minimized_files = multi_minimize(
                self.pdb_files,
                ph=ph,
                force_field=force_field,
                fill_side_chain=fill_side_chain,
                add_missing_loop=add_missing_loop,
                del_water=del_water,
                watdist=watdist,
                rmsd_cutoff=rmsd_cutoff,
                save_dir="minimize",
                overwrite=overwrite,
                cpu_num=self.cpu_num,
            )

        minimized_files = []
        cocrystal_ligands = []
        for minimized_file in self.minimized_files:
            try:
                pdbid = minimized_file.metadata.pdbid
                lig_res_name = minimized_file.metadata.ligand_name
                ligand = minimized_file.find_ligands(specified_names=[lig_res_name])
            except ValueError:
                logger.warning(f"No ligand named {lig_res_name} found in the structure of {pdbid}")
                continue
            if len(ligand) > 1:
                msg = (
                    f"Multiple ligand named {lig_res_name} are found in the structure of {pdbid},\n"
                )
                msg += f"and the first one will be used as the grid center automatically."
                logger.info(msg)

            cocrystal_ligands.append(ligand[0])
            minimized_files.append(minimized_file)
        self.cocrystal_ligands = cocrystal_ligands
        self.minimized_files = minimized_files

        return self.minimized_files

    def _split_cocrystal_lig(self) -> Tuple[list, list]:
        """Split the cocrystal ligand and protein from the minimized structure file.

        Returns:
            Tuple[list, list]: protein_files, ligand_files
        """
        protein_files = []
        ligand_files = []
        with ChDir(self.save_path):
            os.makedirs("protein", exist_ok=True)
            os.makedirs("ligands", exist_ok=True)
            for minimize_file, ligand in zip(self.minimized_files, self.cocrystal_ligands):
                pro_st, lig_st = split_complex(
                    minimize_file, ligand_atom_indexes=ligand.atom_indexes
                )
                protein_file = os.path.join("protein", f"{minimize_file.file_prefix}_pro.pdb")
                ligand_file = os.path.join("ligands", f"{minimize_file.file_prefix}_lig.pdb")
                pro_st.write(protein_file)
                lig_st.write(ligand_file)
                protein_files.append(protein_file)
                ligand_files.append(ligand_file)
        return protein_files, ligand_files

    def grid_generate(
        self,
        box_size: int = 20,
        force_field: Literal["OPLS4", "OPLS3e", "OPLS3", "OPLS_2005"] = "OPLS4",
        overwrite: bool = False,
    ) -> List[GridFile]:
        """Generate grid files for docking.

        Args:
            box_size (int, optional): Box size of the grid. Defaults to 20.
            force_field (Literal['OPLS4', 'OPLS3e', 'OPLS3', 'OPLS_2005'], optional): Force field for the grid. Defaults to 'OPLS4'.
            overwrite (bool, optional): Whether to overwrite existing files. Default to False.

        Raises:
            ValueError: No minimized structure files found.

        Returns:
            List[GridFile]: Grid files list
        """
        if not self.minimized_files:
            raise ValueError("No minimized structure files found.")
        box_size_list = [box_size] * len(self.minimized_files)
        with ChDir(self.save_path):
            self.grid_files = multi_grid_generate(
                self.minimized_files,
                box_center_molnum_list=[lig.molnum for lig in self.cocrystal_ligands],
                box_size_list=box_size_list,
                force_field=force_field,
                save_dir="grid",
                overwrite=overwrite,
                cpu_num=self.cpu_num,
            )
        return self.grid_files

    def dock(
        self,
        retrospective: bool = False,
        force_field: Literal["OPLS4", "OPLS3e", "OPLS3", "OPLS_2005"] = "OPLS4",
        precision: Literal["SP", "XP", "HTVS"] = "SP",
        calc_rmsd: bool = False,
        include_receptor: bool = False,
        overwrite: bool = False,
    ) -> List[dict]:
        """Perform ensemble docking.

        Args:
            retrospective (bool, optional): Whether to add cocrystal molecules to ligands during ensemble docking. Defaults to False.
            force_field (Literal['OPLS4', 'OPLS3e', 'OPLS3', 'OPLS_2005'], optional): Force field for the docking. Defaults to 'OPLS4'.
            precision (Literal['SP', 'XP', 'HTVS'], optional): _description_. Defaults to 'SP'.
            calc_rmsd (bool, optional): Whether to calculate RMSD. Defaults to False.
            include_receptor (bool, optional): Whether to include receptor in the docking. Defaults to False.
            overwrite (bool, optional): Whether to overwrite existing files. Default to False.

        Raises:
            ValueError: No grid files found.

        Returns:
            List[dict]: Docking results
        """
        if not self.grid_files:
            raise ValueError("No grid files found.")

        if retrospective:
            cocrystal_ligand_files = self._split_cocrystal_lig()[1]
            logger.info(
                f"{len(cocrystal_ligand_files)} cocrystal molecules are added to the docking ligands."
            )
            library_files = cocrystal_ligand_files + self.library_files
        docking_pairs = get_docking_pairs(self.grid_files, library_files)
        logger.info(f"{len(docking_pairs)} docking jobs are prepared to run")
        with ChDir(self.save_path):
            self.ensemble_docking_results = multi_dock(
                docking_pairs,
                force_field=force_field,
                precision=precision,
                calc_rmsd=calc_rmsd,
                include_receptor=include_receptor,
                save_dir="dock",
                overwrite=overwrite,
                cpu_num=self.cpu_num,
            )
        self.dock_result_data = multi_extract_data(
            self.ensemble_docking_results, cpu_num=self.cpu_num
        )
        result_file = os.path.join(self.save_path, "dock", f"dock_result_{DATE}.csv")
        self.result_df = pd.DataFrame(self.dock_result_data)
        self.result_df.to_csv(result_file, index=False)
        logger.info(f"Docking result data saved to {result_file}")
        return self.dock_result_data
