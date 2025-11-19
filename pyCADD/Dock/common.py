from typing import List, Union

from pyCADD.utils.common import File

from .const import AMINO_ACIDS, ATOM_RECORDS
from .utils import check_pdb


class PDBLine:
    def __init__(self, line: str) -> None:
        """One line in a pdb file

        Args:
            line (str): one line in a pdb file
        """
        self._line = line.strip()
        self.record_name = ""  # record name
        self.atom_idx = ""  # atom serial number
        self.atom_name = ""  # atom type
        self.alt_loc = ""  # alternate location indicator
        self.res_name = ""  # residue name
        self.chain_id = ""  # chainID
        self.res_id = ""  # resSeq
        self.insertion_code = ""  # iCode
        self.coord_x = ""  # x
        self.coord_y = ""  # y
        self.coord_z = ""  # z
        self.occupancy = ""  # occupancy
        self.temp_factor = ""  # tempFactor
        self.element = ""  # element symbol, right-justified
        self.charge = ""  # charge on the atom

        self._parse()

    @property
    def _slice_define(self):
        return {
            "record_name": (0, 6),
            "atom_idx": (6, 11),
            "atom_name": (12, 16),
            "alt_loc": (16, 17),
            "res_name": (17, 20),
            "chain_id": (21, 22),
            "res_id": (22, 26),
            "insertion_code": (26, 27),
            "coord_x": (30, 38),
            "coord_y": (38, 46),
            "coord_z": (46, 54),
            "occupancy": (54, 60),
            "temp_factor": (60, 66),
            "element": (76, 78),
            "charge": (78, 80),
        }

    def _get_line_slice(self, start: int, end: int, strip: bool = True):
        try:
            info = self._line[start:end]
        except Exception:
            info = ""
        return info.strip() if strip else info

    def _parse(self):
        for key, (start, end) in self._slice_define.items():
            setattr(self, key, self._get_line_slice(start, end))

    def __str__(self) -> str:
        return f"<PDBLine {self.line} >"

    def __repr__(self) -> str:
        return self.__str__()

    @property
    def is_atom_line(self) -> bool:
        return self._line.startswith("ATOM") or self._line.startswith("HETATM")

    @property
    def is_amino(self) -> bool:
        return self.record_name in ATOM_RECORDS and self.res_name in AMINO_ACIDS

    @property
    def is_hetatm(self) -> bool:
        return self._line.startswith("HETATM")

    @property
    def is_conect(self) -> bool:
        return self._line.startswith("CONECT")

    @property
    def is_ter(self) -> bool:
        return self._line.startswith("TER")

    @property
    def is_end(self) -> bool:
        return self._line.startswith("END")

    @property
    def _formatter(self) -> str:
        if len(self.get_atom_name()) == 4:
            return "{:<6}{:>5} {:<4}{:1}{:<3} {:1}{:>4}{:1}   {:>8}{:>8}{:>8}{:>6}{:>6}          {:>2}{:>2}"
        # start writing atom name at 14th column
        return "{:<6}{:>5}  {:<3}{:1}{:<3} {:1}{:>4}{:1}   {:>8}{:>8}{:>8}{:>6}{:>6}          {:>2}{:>2}"

    @property
    def line(self) -> str:
        return self._formatter.format(
            self.record_name,
            self.atom_idx,
            self.atom_name,
            self.alt_loc,
            self.res_name,
            self.chain_id,
            self.res_id,
            self.insertion_code,
            self.coord_x,
            self.coord_y,
            self.coord_z,
            self.occupancy,
            self.temp_factor,
            self.element,
            self.charge,
        )

    def get_line(self) -> str:
        """Get the line string from current attributes

        Returns:
            str: line string
        """
        return self.line

    def get_atom_name(self) -> str:
        """Get the atom name from the line

        Returns:
            str: atom name
        """
        atom_name = self.atom_name
        if not atom_name:
            return atom_name
        elif atom_name[0].isdigit():  # for name such as: 1HB, 1HG2
            atom_name = atom_name[1:] + atom_name[0]
        return atom_name


class PDBLineParser:
    def __init__(self, pdb_str: str = None, pdb_file: str = None) -> None:
        """Parse pdb file string or file.

        Args:
            pdb_str (str, optional): pdb string. Defaults to None.
            pdb_file (str, optional): pdb file path. Required if pdb_str is not provided. Defaults to None.

        Raises:
            ValueError: Either pdb_str or pdb_file must be provided.
        """
        if pdb_str is None and pdb_file is None:
            raise ValueError("Either pdb_str or pdb_file must be provided")
        self.pdb_file = pdb_file
        self.pdb_str = pdb_str if pdb_str is not None else self._read_pdb_file()
        self.pdb_lines = []
        self._idx_map = {}
        self._parse_lines()

    def __str__(self) -> str:
        return "\n".join(self.get_line_str_list())

    def __repr__(self) -> str:
        return self.__str__()

    def _parse_lines(self):
        self.pdb_lines = [PDBLine(line) for line in self.pdb_str.splitlines() if line.strip()]
        self._idx_map = {line.atom_idx: i for i, line in enumerate(self.pdb_lines)}

    def _read_pdb_file(self):
        """Read the pdb file if pdb_str is not provided

        Returns:
            str: pdb file content
        """
        with open(self.pdb_file) as f:
            _pdb_str = f.read()
        return _pdb_str

    def get_lines(self):
        """Get a list of all lines parsed from the pdb file.

        Returns:
            list[PDBLine]: pdb line object list
        """
        return self.pdb_lines

    def get_atom_lines(self):
        """Get a list of atom line objects parsed from the pdb file.

        Returns:
            list[PDBLine]: atom line object list
        """
        return [line for line in self.pdb_lines if line.is_atom_line]

    def get_amino_lines(self):
        """Get a list of amino acid line objects parsed from the pdb file.
        Amino acid line is defined as the line with record name 'ATOM' and res_name in AMINO_ACIDS

        Returns:
            list[PDBLine]: amino acid line objects list
        """
        return [line for line in self.pdb_lines if line.is_amino]

    def get_hetatm_lines(self):
        """Get a list of HETATM line objects parsed from the pdb file.

        Returns:
            list[PDBLine]: HETATM line objects list
        """
        return [line for line in self.pdb_lines if line.is_hetatm]

    def get_str_list(self):
        """Get a list of pdb line strings.

        Returns:
            list[str]: pdb line strings list
        """
        return [line.get_line() for line in self.pdb_lines]

    def save_pdb(self, file_path):
        """Save the pdb file to the file_path

        Args:
            file_path (str): file path to save the pdb file
        """
        with open(file_path, "w") as f:
            f.write("\n".join(self.get_line_str_list()))


class PDBFile(File):
    def __init__(self, path: str) -> None:
        """PDB file class

        Args:
            path (str): file path string
        """
        super().__init__(path)
        self.pdbid = self.file_prefix[:4] if check_pdb(self.file_prefix[:4]) else None
        self.pdb_parser = PDBLineParser(pdb_file=self.file_path)

    def _catch_lig(self) -> list:
        result_list = []
        _items = ["id", "chain", "resid", "atom_num"]
        with open(self.file_path, "r") as f:
            lines = f.read().splitlines()
        for line in lines:
            if line.startswith("HET "):
                match = ",".join(line.split()[1:])
                lig_dict = {k: v for k, v in zip(_items, match.split(","))}
                result_list.append(lig_dict)
        return result_list

    def get_lines(self, return_str: bool = False) -> Union[list, str]:
        """Get the pdb file content as a list of lines or a string

        Args:
            return_str (bool, optional): get the string instead of list[str]. Defaults to False.

        Returns:
            list|str: pdb file content
        """
        return (
            self.pdb_parser.get_str_list()
            if not return_str
            else "\n".join(self.pdb_parser.get_str_list())
        )

    def get_chain(self, chain_id: str, return_str: bool = False) -> Union[list, str]:
        """Get the pdb file content of a single chain

        Args:
            chain_id (str): chain id
            return_str (bool, optional): get the string instead of list[str]. Defaults to False.

        Returns:
            str: pdb file content of a single chain
        """
        chain_content = [
            line.get_line()
            for line in self.pdb_parser.get_atom_lines()
            if line.chain_id == chain_id
        ]
        return "\n".join(chain_content) if return_str else chain_content


class EnsembleInputFile(File):
    def __init__(self, path: str) -> None:
        """Input file class for ensemble docking

        Args:
            path (_type_): file path
        """
        super().__init__(path)
        self.mappings = None
        self._pdbid_list = None
        self._ligand_list = None

    @property
    def pdbid_list(self) -> list:
        if self._pdbid_list is None:
            self._pdbid_list = self.get_pdbid_list()
        return self._pdbid_list

    @property
    def ligand_list(self) -> list:
        if self._ligand_list is None:
            self._ligand_list = self.get_ligand_list()
        return self._ligand_list

    @classmethod
    def from_csv(cls, file_path: str, sep: str = ",", header: bool = False) -> "EnsembleInputFile":
        """
        Parse input file as csv format

        Args:
            file_path (str): csv file path
            sep (str, optional): separator. Defaults to ','.
            header (bool, optional): whether the csv file has header. Defaults to None.

        csv examples:
            ```
            1XJ7,DHT
            1XQ3,R18
            2AM9,TES
            2AM9,DTT
            2YLP,TES
            2YLP,056
            ```

        Returns:
            EnsembleInputFile: instance of EnsembleInputFile
        """
        csv_file = File(file_path)
        with open(csv_file.file_path, "r") as f:
            raw_list = f.read().splitlines()
        if header:
            raw_list = raw_list[1:]

        mappings = []
        for line in raw_list:
            item = line.split(sep)
            if len(item) == 1:
                pdbid = line[0].strip()
                ligand_name = ""
            elif len(line) >= 2:
                pdbid, ligand_name = item[0].strip(), item[1].strip()

            mappings.append({"receptor": csv_file.file_prefix, "pdb": pdbid, "ligand": ligand_name})
        ins = cls(file_path)
        ins.mappings = mappings
        return ins

    @classmethod
    def from_ini(cls, config_file: str) -> "EnsembleInputFile":
        """Parse input file as ini format

        Args:
            config_file (str): ini file path

        ini examples:
            ```
            [P10275]
                1XJ7: DHT
                1XQ3: R18
                2AM9: TES,DTT
                2YLP: TES,056
            ```

        Returns:
            EnsembleInputFile: instance of EnsembleInputFile
        """
        from pyCADD.utils.common import FixedConfig

        config = FixedConfig()
        config.read(config_file)
        receptors = [receptor for receptor in config.sections()]
        mappings = []
        for receptor in receptors:
            for _item in config.items(receptor):
                ligs = _item[1].split(",")
                for lig in ligs:
                    mappings.append({"receptor": receptor, "pdb": _item[0], "ligand": lig})
        ins = cls(config_file)
        ins.mappings = mappings
        return ins

    @classmethod
    def from_yaml(cls, yaml_file: str) -> "EnsembleInputFile":
        """Parse input file as yaml format

        Args:
            yaml_file (str): yaml file path


        yaml examples:
            ```
            P10275:
                1XJ7: DHT
                1XQ3:
                - R18
                2AM9:
                - TES
                - DTT
                2YLP:
                - TES
                - '056'

        Returns:
            EnsembleInputFile: instance of EnsembleInputFile
        """
        import yaml

        with open(yaml_file, "r") as f:
            yaml_dict = yaml.load(f, Loader=yaml.FullLoader)

        mappings = []
        for receptor in yaml_dict.keys():
            for pdb, ligs in yaml_dict[receptor].items():
                if isinstance(ligs, str):
                    ligs = [ligs]
                for lig in ligs:
                    mappings.append({"receptor": receptor, "pdb": pdb, "ligand": lig})
        ins = cls(yaml_file)
        ins.mappings = mappings
        return ins

    @classmethod
    def parse_file(cls, path: str, header: bool = False) -> "EnsembleInputFile":
        """Parse input file

        Args:
            path (str): file path
            header (bool, optional): whether the file has header. Only for csv file. Defaults to False.

        Raises:
            ValueError: Unsupported file type

        Returns:
            EnsembleInputFile: instance of EnsembleInputFile
        """
        file = File(path)
        if file.file_ext.lower() in ["csv", "txt"]:
            return cls.from_csv(path, header=header)
        elif file.file_ext.lower() == "ini":
            return cls.from_ini(path)
        elif file.file_ext.lower() in ["yaml", "yml"]:
            return cls.from_yaml(path)
        else:
            raise ValueError(f"Unsupported file type: {file.file_path}")

    def read(self, file_path: str) -> None:
        """Read and parse the input file

        Args:
            file_path (str): file path
        """
        self.mappings = self.parse_file(file_path).mappings

    def get_pairs_list(self) -> List[tuple]:
        """Get the list of pairs. Pairs are defined as (pdb, ligand)

        Returns:
            list: list of pairs
        """
        if self.mappings is None:
            self.read(self.file_path)
        return [(item["pdb"], item["ligand"]) for item in self.mappings]

    def get_pdbid_list(self) -> List[str]:
        """Get the list of unique pdb ids

        Returns:
            list: list of pdb ids
        """
        if self.mappings is None:
            self.read(self.file_path)
        return sorted(set([item["pdb"] for item in self.mappings]))

    def get_ligand_list(self) -> List[str]:
        """Get the list of unique ligands

        Returns:
            list: list of ligands
        """
        if self.mappings is None:
            self.read(self.file_path)
        return sorted(set([item["ligand"] for item in self.mappings]))
