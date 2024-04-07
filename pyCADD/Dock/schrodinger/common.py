import os
from dataclasses import dataclass
from typing import List, Union

from pyCADD.utils.common import File

from . import Structure, StructureReader, struc


@dataclass
class MetaData:
    """Data class for metadata"""
    pdbid: str = 'Undefined'
    ligand_name: str = 'Undefined'
    ligand_resnum: int = 0
    _action: str = 'Undefined'
    docking_ligand_name: str = 'Undefined'
    precision: str = 'Undefined'
    
    property_index_map = {
        "pdbid": 0,
        "ligand_name": 1,
        "_action": 2,
        "docking_ligand_name": 3,
        "precision": 4
    }

    @property
    def internal_ligand_name(self) -> str:
        return self.ligand_name

    @property
    def internal_ligand_resnum(self) -> int:
        return self.ligand_resnum

    @classmethod
    def parse_from_filename(cls, file_name: str, sep: str = '_') -> 'MetaData':
        """Parse metadata from the filename

        Args:
            file_name (str): file name or path
            sep (str, optional): separator. Defaults to '_'.
            ligand_resnum (int, optional): ligand residue number. Defaults to 0.
            kwargs: other attributes to be added

        Returns:
            MetaData: metadata object
        """
        file_name, _ = os.path.splitext(os.path.basename(file_name))
        metadatas = file_name.split(sep=sep)
        ins = cls()
        for p, index in cls.property_index_map.items():
            if index < len(metadatas):
                setattr(ins, p, metadatas[index])
        return ins

    def generate_file_name(self, sep: str = "_") -> str:
        """Generate file name from metadata

        Returns:
            str: file name without extension
        """
        return f"{self.pdbid}{sep}{self.ligand_name}"
    
    def __str__(self) -> str:
        return f"{self.__dict__}"
    
    def __repr__(self) -> str:
        return self.__str__()


class MaestroFile(File):
    def __init__(self, path: str, **kwargs) -> None:
        """Mastro file class

        Args:
            path (str): file path
        """
        super().__init__(path)
        self.metadata = MetaData.parse_from_filename(self.file_path, **kwargs)

    @property
    def st_reader(self) -> StructureReader:
        return StructureReader(self.file_path)

    @property
    def structures(self) -> List[Structure]:
        return [st for st in self.st_reader]

    def __str__(self) -> str:
        return f"<Maestro File at {self.file_path} with {len(self.structures)} structure(s)>\n{self.metadata}"
    
    def __repr__(self) -> str:
        return self.__str__()
    
    @staticmethod
    def get_structure(file_path: str, index: int = 0) -> Structure:
        """Get the structure from the file

        Args:
            file_path (str): maestro file path
            index (int, optional): index of the structure to get. Defaults to 0.

        Returns:
            Structure: Maestro structure object
        """
        return MaestroFile(file_path).structures[index]

    def get_residue(self, resnum: int, structure_index: int = 0) -> Union[struc._Residue, None]:
        """Get the Residue object from the structure

        Args:
            mol_resnum (int): qeury molecule residue number
            structure_index (int): index of the structure

        Returns:
            struc.Residue: Maestro Residue object
        """
        return next((res for res in self.structures[structure_index].residues if res.resnum == resnum), None)

    def get_molecule_by_res(self, mol_resnum: int = None, structure_index: int = 0) -> Union[struc._Molecule, None]:
        """Get the Molecule object from the structure by residue number.

        Args:
            mol_resnum (int): qeury molecule residue number
            structure_index (int): index of the structure

        Returns:
            struc._Molecule: Maestro Molecule object
        """
        res = self.get_residue(mol_resnum, structure_index)
        if res is None:
            return None
        return self.structures[structure_index].molecules[res.molecule_number]

    def get_covalent_bond(self, mol_resnum: int, structure_index: int = 0) -> Union[list, None]:
        """
        Get list of covalent bond(s) between query molecule and other residues. 
        Return None if no covalent bond is found.

        Args:
            mol_resnum (int): query residue number of the molecule
            structure_index (int): index of the structure
        """
        bonds = self.structures[structure_index].bond
        covalent_bonds = []
        for bond in bonds:
            resnum1 = bond.atom1.getResidue().resnum
            resnum2 = bond.atom2.getResidue().resnum
            if resnum1 != resnum2 and mol_resnum in (resnum1, resnum2):
                covalent_bonds.append(bond)

        if len(covalent_bonds) == 0:
            return None
        return covalent_bonds

    def get_molnum_by_res(self, mol_resnum: int, structure_index: int = 0) -> int:
        """Get the molecule number from the structure by residue number.

        Args:
            mol_resnum (int, optional): query residue number. Defaults to None.
            structure_index (int): index of the structure

        Returns:
            int: Molecule number
        """
        res = self.get_residue(mol_resnum, structure_index)
        if res is None:
            raise ValueError(
                f"Residue {mol_resnum} not found in the structure.")
        return res.molecule_number


class DockResultFile(MaestroFile):
    def __init__(self, path: str, ligand_only: bool) -> None:
        """Docking result file class

        Args:
            path (str): file path
            ligand_only (bool): True if the file contains only ligand structure
        """
        super().__init__(path)
        self.ligand_only = ligand_only

    def get_receptor_structure(self) -> Structure:
        """Get the receptor structure from the docking result file

        Returns:
            Structure: Maestro structure object
        """
        if self.ligand_only or len(self.structures == 1):
            return None
        return self.structures[0]

    def get_ligand_structure(self) -> Structure:
        """Get the docking ligand structure from the docking result file

        Returns:
            Structure: Maestro structure object
        """
        if self.ligand_only:
            return self.structures[0]
        return self.structures[1]

    def get_result_dict(self) -> dict:
        """Get the docking result information

        Returns:
            dict: docking result information
        """
        return {k: v for k, v in self.get_ligand_structure().property.item()}
