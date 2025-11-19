import os
from copy import deepcopy
from dataclasses import dataclass
from typing import List, Union, Any

from pyCADD.utils.common import File

from . import AslLigandSearcher, Ligand, Structure, StructureReader, struc
from .config import DataConfig

DEBUG = os.getenv("PYCADD_DEBUG", False)


@dataclass
class MetaData:
    """Data class for metadata"""

    pdbid: str = ""
    ligand_name: str = ""
    action: str = ""
    docking_ligand_name: str = ""
    precision: str = ""
    _PARSE_DICT = {
        "pdbid": 0,
        "ligand_name": 1,
        "action": 2,
        "docking_ligand_name": 3,
        "precision": 4,
    }

    @property
    def internal_ligand_name(self) -> str:
        return self.ligand_name

    @classmethod
    def parse_from_filename(
        cls, file_name: str, sep: str = "_", parse_dict: dict = None
    ) -> "MetaData":
        """Parse metadata from the filename

        Args:
            file_name (str): file name or path
            sep (str, optional): separator. Defaults to '_'.
            parse_dict (dict, optional): specify position index that metadata parse attrs from the file name. Defaults to {'pdbid': 0, 'ligand_name': 1, 'action': 2, 'docking_ligand_name': 3, 'precision': 4}.

        Returns:
            MetaData: metadata object
        """
        position_index_map = parse_dict if parse_dict is not None else cls._PARSE_DICT
        file_name, _ = os.path.splitext(os.path.basename(file_name))
        metadatas = file_name.split(sep=sep)
        ins = cls()
        for attr, pos in position_index_map.items():
            if pos < len(metadatas):
                setattr(ins, attr, metadatas[pos])
        return ins

    @classmethod
    def parse_from_dict(cls, metadata_dict: dict) -> "MetaData":
        """Parse metadata from the dictionary

        Args:
            metadata_dict (dict): metadata dictionary

        Returns:
            MetaData: metadata object
        """
        ins = cls()
        for key, value in metadata_dict.items():
            if hasattr(ins, key):
                setattr(ins, key, value)
        return ins

    @classmethod
    def parse_from_metadata(cls, metadata: "MetaData") -> "MetaData":
        """Parse metadata from the metadata object

        Args:
            metadata (MetaData): metadata object

        Returns:
            MetaData: metadata object
        """
        return metadata.copy()

    def generate_file_name(
        self,
        attributes: list = ["pdbid", "ligand_name", "action", "docking_ligand_name", "precision"],
        sep="_",
    ) -> str:
        """Generate file name from metadata

        Args:
            attributes (list, optional): attributes to be included in the file name. Defaults to ['pdbid', 'ligand_name', 'action', 'docking_ligand_name', 'precision'].
            sep (str, optional): separator. Defaults to '_'.

        Returns:
            str: file name without extension

        Example:
            >>> metadata = MetaData(pdbid='1ABC', ligand_name='LIG', action='glide-dock', docking_ligand_name='LIG', precision='SP')
            >>> metadata.generate_file_name()
            '1ABC_LIG_glide-dock_LIG_SP'
            >>> metadata.generate_file_name(attributes=['pdbid', 'ligand_name', 'action'], sep='-')
            '1ABC-LIG-glide-dock'
        """
        for attr in attributes:
            if not hasattr(self, attr):
                raise ValueError(f"Attribute {attr} is not found in the metadata.")
        return sep.join([str(getattr(self, a)) for a in attributes])

    def __str__(self) -> str:
        return f"{self.__dict__}"

    def __repr__(self) -> str:
        return self.__str__()

    def copy(self) -> "MetaData":
        """Get the deep copy of the metadata

        Returns:
            MetaData: metadata object
        """
        return deepcopy(self)

    def set(self, attr, value) -> None:
        """Set the attribute value

        Args:
            attr (str): attribute name
            value (any): attribute value
        """
        setattr(self, attr, value)

    class NoDefault:
        pass

    def get(self, attr, default=NoDefault) -> Any:
        """Get the attribute value

        Args:
            attr (str): attribute name
            default (any): default value that will be returned if the attribute is not found; without it, an AttributeError will be raised in this case.

        Returns:
            any: attribute value. If the attribute is not found, return provided default value.
        """
        if default is self.NoDefault:
            return getattr(self, attr)
        return getattr(self, attr, default)

    def delete(self, attr) -> None:
        """Delete the attribute

        Args:
            attr (str): attribute name
        """
        delattr(self, attr)


class BaseMaestroFile(File):
    def __init__(self, path: str, metadata: MetaData = None, **kwargs) -> None:
        """Base Mastro file class

        Args:
            path (str): file path
            metadata (MetaData, optional): metadata object. Defaults to None.
        """
        if isinstance(path, BaseMaestroFile):
            super().__init__(path.file_path, **kwargs)
            self.metadata = path.metadata.copy()
        else:
            super().__init__(path, **kwargs)
            if metadata is None:
                self.metadata = MetaData.parse_from_filename(self.file_path)
            else:
                self.metadata = MetaData.parse_from_metadata(metadata)


class GridFile(BaseMaestroFile):
    def __init__(self, path: str, metadata: MetaData = None, **kwargs) -> None:
        """Mastro grid file class

        Args:
            path (str): file path
            metadata (MetaData, optional): metadata object. Defaults to None.
        """
        super().__init__(path, metadata, **kwargs)


class LigandSearched(Ligand):
    def __init__(self, ligand_obj: Ligand):
        self.__dict__ = ligand_obj.__dict__
        self._ligand = ligand_obj

    @property
    def pdbres(self) -> str:
        return self._ligand.pdbres.strip().upper()

    @property
    def chain(self) -> str:
        return self._ligand.atom_objects[0].chain

    @property
    def resnum(self) -> int:
        return self._ligand.atom_objects[0].getResidue().resnum

    @property
    def molnum(self) -> int:
        return self.mol_num

    def __str__(self):
        return f"<LigandSearched: name={self.pdbres.strip()} chain={self.chain} resnum={self.resnum} molnum={self.molnum}>"

    def __repr__(self):
        return self.__str__()


class MaestroFile(BaseMaestroFile):

    def __init__(self, path: str, metadata: MetaData = None, **kwargs) -> None:
        """Maestro file class

        Args:
            path (str): file path
            metadata (MetaData, optional): metadata object. Defaults to None.
        """
        super().__init__(path, metadata, **kwargs)
        self._ligands = []

    @property
    def st_reader(self) -> StructureReader:
        if not self.file_path or not os.path.exists(self.file_path):
            return []
        return StructureReader(self.file_path)

    @property
    def structures(self) -> List[Structure]:
        return [st for st in self.st_reader]

    def __str__(self) -> str:
        info = f"<Maestro File at {self.file_path} with {len(self.structures)} structure(s)>"
        if DEBUG:
            info += f"\n{self.metadata}"
        return info

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

    def get_chain_structure(self, chain_id: str, structure_index: int = 0) -> Structure:
        """Get the chain structure from the structure file

        Args:
            chain_id (str): chain ID
            structure_index (int): index of the structure

        Returns:
            Structure: Maestro structure object
        """
        return self.structures[structure_index].chain[chain_id].extractStructure()

    def get_residue(self, resnum: int, structure_index: int = 0) -> Union[struc._Residue, None]:
        """Get the Residue object from the structure

        Args:
            mol_resnum (int): qeury molecule residue number
            structure_index (int): index of the structure

        Returns:
            struc.Residue: Maestro Residue object
        """
        return next(
            (res for res in self.structures[structure_index].residue if res.resnum == resnum), None
        )

    def get_molecule_by_res(
        self, resnum: int = None, structure_index: int = 0
    ) -> Union[struc._Molecule, None]:
        """Get the Molecule object from the structure by residue number.

        Args:
            resnum (int): qeury molecule residue number
            structure_index (int): index of the structure

        Returns:
            struc._Molecule: Maestro Molecule object
        """
        res = self.get_residue(resnum, structure_index)
        if res is None:
            return None
        return self.structures[structure_index].molecule[res.molecule_number]

    def get_covalent_bond(self, resnum: int, structure_index: int = 0) -> Union[list, None]:
        """
        Get list of covalent bond(s) between query molecule and other residues.
        Return None if no covalent bond is found.

        Args:
            resnum (int): query residue number of the molecule
            structure_index (int): index of the structure

        Returns:
            list: list of covalent bond(s) between query molecule and other residues
        """
        bonds = self.structures[structure_index].bond
        covalent_bonds = []
        for bond in bonds:
            resnum1 = bond.atom1.getResidue().resnum
            resnum2 = bond.atom2.getResidue().resnum
            if resnum1 != resnum2 and resnum in (resnum1, resnum2):
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
            raise ValueError(f"Residue {mol_resnum} not found in the structure.")
        return res.molecule_number

    def find_ligands(
        self,
        specified_names: list = None,
        included_names: list = None,
        excluded_names: list = None,
        min_heavy_atom_count: int = None,
        max_atom_count: int = None,
        allow_amino_acid_only_molecules: bool = False,
        allow_ion_only_molecules: bool = False,
        structure_index: int = 0,
    ) -> List[LigandSearched]:
        """Find all ligands from the structure file

        Args:
            specified_names (list, optional): PDB residue names which will be considered as ligands. \
                Any found residues not in this list will be filtered out if this list is not None. Defaults to None.
            included_names (list, optional): PDB residue names which always be considered as ligands. Defaults to None.
            excluded_names (list, optional): PDB residue names which will never be considered as ligands. Defaults to None.
            min_heavy_atom_count (int, optional): Minimum number of heavy atoms required in each ligand molecule. Defaults to None.
            max_atom_count (int, optional): Maximum number of heavy atoms for a ligand molecule (does not include hydrogens). Defaults to None.
            allow_amino_acid_only_molecules (bool, optional): If True, consider small molecules containing only amino acids to be ligands. Defaults to False.
            allow_ion_only_molecules (bool, optional): If True, Consider charged molecules to be ligands. Defaults to False.
            structure_index (int): Index of the structure

        Returns:
            list[LigandSearched]: List of found Ligand objects.
                useful properties: pdbres, chain, resnum, mol_num, centriod, st, atom_indexes
        """
        included_names = set(included_names) if included_names is not None else set()
        excluded_names = set(excluded_names) if excluded_names is not None else set()
        allow_amino_acid_only_molecules = bool(allow_amino_acid_only_molecules)
        allow_ion_only_molecules = bool(allow_ion_only_molecules)
        args = {
            "included_residue_names": included_names,
            "excluded_residue_names": excluded_names,
            "allow_amino_acid_only_molecules": allow_amino_acid_only_molecules,
            "allow_ion_only_molecules": allow_ion_only_molecules,
        }
        if min_heavy_atom_count is not None:
            args["min_heavy_atom_count"] = int(min_heavy_atom_count)
        if max_atom_count is not None:
            args["max_atom_count"] = int(max_atom_count)

        searcher = AslLigandSearcher(**args)
        self._ligands = searcher.search(self.structures[structure_index])
        self._ligands = [LigandSearched(lig) for lig in self._ligands]
        if specified_names is not None:
            self._ligands = [lig for lig in self._ligands if lig.pdbres in specified_names]
            if not self._ligands:
                raise ValueError(f"No ligand named {specified_names} found in {self.file_path}")
        return self._ligands


class DockResultFile(MaestroFile):
    def __init__(
        self, path: str, metadata: MetaData = None, include_receptor: bool = None, **kwargs
    ) -> None:
        """Docking result file class

        Args:
            path (str): file path
            metadata (MetaData, optional): metadata object. Defaults to None.
            include_receptor (bool): True if the file contains receptor structure.
        """
        super().__init__(path, **kwargs)
        self.include_receptor = False
        if isinstance(metadata, MetaData):
            self.metadata = metadata.copy()
            self.include_receptor = metadata.get("include_receptor", False)
        self.include_receptor = (
            include_receptor if include_receptor is not None else self.include_receptor
        )

    def get_receptor_structure(self) -> Structure:
        """Get the receptor structure from the docking result file

        Returns:
            Structure: Maestro structure object
        """
        if not self.include_receptor:
            return None
        return self.structures[0]

    def get_ligand_structure(self) -> Structure:
        """Get the docking ligand structure from the docking result file

        Returns:
            Structure: Maestro structure object
        """
        if self.include_receptor:
            return self.structures[1]
        return self.structures[0]

    def get_raw_results(self) -> List[dict]:
        """Get the raw docking result information of all structures

        Returns:
            list[dict]: raw docking result information
        """
        return [{k: v for k, v in st.property.items()} for st in self.structures]

    def get_results(self) -> List[dict]:
        """Get the docking result information of all structures

        Returns:
            list[dict]: docking result information
        """
        config = DataConfig(precision=self.metadata.precision)
        result_dict = {
            "pdbid": self.metadata.pdbid,
            "precision": self.metadata.precision,
            "internal_ligand_name": self.metadata.internal_ligand_name,
            "docking_ligand_name": self.metadata.docking_ligand_name,
        }
        results = []
        for st in self.structures:
            _result_dict = result_dict.copy()
            _result_dict.update(
                {key: st.property.get(value, None) for key, value in config.properties.items()}
            )
            results.append(_result_dict)
        return results
