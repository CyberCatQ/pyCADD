"""Core utilities for Gaussian quantum chemistry calculations and molecular density analysis.

This module provides core functions for generating Gaussian input files, processing molecular
structures, and handling charge calculations for quantum chemistry workflows.
"""

import logging
import os
import tempfile

from pyCADD.utils.common import ChDir, File
from pyCADD.utils.tool import read_file, shell_run

logger = logging.getLogger(__name__)

# MOL2 file format string for atom line formatting
MOL2_LINE_FMT = "{:>7} {:<10} {:>10} {:>10} {:>10} {:<5} {:>4} {:<7} {:>10}"

# Template script for generating Gaussian input files using OpenBabel
GAUSS_SCRIPT_TEMP = """
cat << EOF > {gauss_file}
%chk={chk_file}
%mem={mem_use}
%nproc={cpu_num}
{keyword}

{title}

  {charge} {multiplicity}
EOF
obabel -i{file_ext} {mol_file} -ogjf | awk '{{if (NR>5 && $1 !~ "[0-9]") print }}' >> {gauss_file}
"""


def generate_gauss_input(
    structure_file: str | File,
    keyword: str,
    chk_file: str = None,
    title: str = None,
    charge: int = 0,
    multiplicity: int = 1,
    mem_use: str = "4GB",
    cpu_num: int = 4,
):
    """Generates a Gaussian input file from a molecular structure.

    Creates a Gaussian input (.gjf) file by combining molecular geometry from
    OpenBabel with specified quantum chemistry parameters.

    Args:
        structure_file (str | File): Path to the input molecular structure file or File object.
        keyword (str): Gaussian calculation keywords (e.g., "# opt B3LYP/6-31G*").
        chk_file (str, optional): Name of the checkpoint file. Defaults to structure_file prefix + ".chk".
        title (str, optional): Title for the calculation. Defaults to structure_file prefix + " Gaussian Input".
        charge (int, optional): Molecular charge. Defaults to 0.
        multiplicity (int, optional): Spin multiplicity. Defaults to 1.
        mem_use (str, optional): Memory allocation string. Defaults to "4GB".
        cpu_num (int, optional): Number of processors to use. Defaults to 4.

    Returns:
        str: Content of the generated Gaussian input file.

    Raises:
        FileNotFoundError: If the structure file doesn't exist.
        RuntimeError: If OpenBabel conversion fails.
    """
    file = File(structure_file)
    chk_file = chk_file or f"{file.file_prefix}.chk"
    title = title or f"{file.file_prefix} Gaussian Input"
    with tempfile.TemporaryDirectory() as tmpdir:
        with ChDir(tmpdir):
            shell_run(
                GAUSS_SCRIPT_TEMP.format(
                    file_ext=file.file_ext.replace("out", "g16"),
                    mol_file=file.file_path,
                    gauss_file="tmp.gjf",
                    chk_file=chk_file,
                    keyword=keyword,
                    title=title,
                    charge=charge,
                    multiplicity=multiplicity,
                    mem_use=mem_use,
                    cpu_num=cpu_num,
                )
            )
            return read_file("tmp.gjf")


def generate_opt_input(
    structure_file: str | File,
    charge: int,
    multiplicity: int = 1,
    dft: str = "B3LYP",
    basis_set: str = "6-31g*",
    solvent: str = "water",
    loose: bool = True,
    dispersion_correct: bool = True,
    td: bool = False,
    freq: bool = False,
    mem_use: str = "4GB",
    cpu_num: int = 4,
):
    """Generates Gaussian input for geometry optimization calculations.

    Creates a Gaussian input file optimized for geometry optimization with
    customizable DFT method, basis set, and calculation options.

    Args:
        structure_file (str | File): Path to the input molecular structure file or File object.
        charge (int): Molecular charge (required parameter).
        multiplicity (int, optional): Spin multiplicity. Defaults to 1.
        dft (str, optional): DFT functional name. Defaults to "B3LYP".
        basis_set (str, optional): Basis set specification. Defaults to "6-31g*".
        solvent (str, optional): Solvent for implicit solvation model. Defaults to "water".
        loose (bool, optional): Whether to use loose convergence criteria. Defaults to True.
        dispersion_correct (bool, optional): Whether to apply dispersion correction (GD3BJ). Defaults to True.
        td (bool, optional): Whether to perform time-dependent DFT calculation. Defaults to False.
        freq (bool, optional): Whether to calculate vibrational frequencies. Defaults to False.
        mem_use (str, optional): Memory allocation string. Defaults to "4GB".
        cpu_num (int, optional): Number of processors to use. Defaults to 4.

    Returns:
        str: Content of the generated Gaussian optimization input file.

    Note:
        The checkpoint file is automatically named as "{structure_prefix}_opt.chk" or
        "{structure_prefix}_opt_td.chk" if TD-DFT is enabled.
    """
    mol_file = File(structure_file)
    mol_name = mol_file.file_prefix
    td_suffix = "_td" if td else ""
    chk_file_path = f"{mol_name}_opt{td_suffix}.chk"

    scrf_kw = f" scrf(solvent={solvent})" if solvent else ""
    opt_kw = "opt=loose" if loose else "opt"
    td_kw = " TD" if td else ""
    correct_kw = " em=GD3BJ" if dispersion_correct else ""
    freq_kw = " freq" if freq else ""
    keyword = f"# {opt_kw} {dft}/{basis_set}{correct_kw}{td_kw}{scrf_kw}{freq_kw}"
    title = f"{mol_name}{td_suffix} Optimization"
    return generate_gauss_input(
        structure_file=structure_file,
        keyword=keyword,
        chk_file=chk_file_path,
        mem_use=mem_use,
        cpu_num=cpu_num,
        title=title,
        charge=charge,
        multiplicity=multiplicity,
    )


def generate_energy_input(
    structure_file: str | File,
    charge: int,
    multiplicity: int,
    dft: str = "B3LYP",
    basis_set: str = "6-31g*",
    solvent: str = "water",
    dispersion_correct: bool = True,
    td: bool = False,
    esp_calculate: bool = False,
    mem_use: str = "4GB",
    cpu_num: int = 4,
):
    """Generates Gaussian input for single-point energy calculations.

    Creates a Gaussian input file for single-point energy calculations with
    options for electrostatic potential (ESP) calculation and TD-DFT.

    Args:
        structure_file (str | File): Path to the input molecular structure file or File object.
        charge (int): Molecular charge (required parameter).
        multiplicity (int): Spin multiplicity (required parameter).
        dft (str, optional): DFT functional name. Defaults to "B3LYP".
        basis_set (str, optional): Basis set specification. Defaults to "6-31g*".
        solvent (str, optional): Solvent for implicit solvation model. Defaults to "water".
        dispersion_correct (bool, optional): Whether to apply dispersion correction (GD3BJ). Defaults to True.
        td (bool, optional): Whether to perform time-dependent DFT calculation. Defaults to False.
        esp_calculate (bool, optional): Whether to calculate electrostatic potential for RESP charges. Defaults to False.
        mem_use (str, optional): Memory allocation string. Defaults to "4GB".
        cpu_num (int, optional): Number of processors to use. Defaults to 4.

    Returns:
        str: Content of the generated Gaussian single-point energy input file.

    Note:
        When esp_calculate is True, the calculation includes "pop=MK IOp(6/33=2,6/42=6)"
        keywords for Merz-Kollman population analysis suitable for RESP charge fitting.
    """
    mol_file = File(structure_file)
    mol_name = mol_file.file_prefix
    td_suffix = "_td" if td else ""
    chk_file_path = f"{mol_name}_energy{td_suffix}.chk"
    scrf_kw = f" scrf(solvent={solvent})" if solvent else ""
    td_kw = " TD" if td else ""
    correct_kw = " em=GD3BJ" if dispersion_correct else ""
    esp_kw = " pop=MK IOp(6/33=2,6/42=6)" if esp_calculate else ""
    keyword = f"# {dft}/{basis_set}{correct_kw}{td_kw}{scrf_kw}{esp_kw}"
    title = f"{mol_name}{td_suffix} Single Point Energy"
    return generate_gauss_input(
        structure_file=structure_file,
        keyword=keyword,
        chk_file=chk_file_path,
        mem_use=mem_use,
        cpu_num=cpu_num,
        title=title,
        charge=charge,
        multiplicity=multiplicity,
    )


def generate_fchk(chk_file_path: str | File):
    """Converts a Gaussian checkpoint file to formatted checkpoint file.

    Uses the Gaussian formchk utility to convert a binary checkpoint (.chk) file
    to a formatted checkpoint (.fchk) file that can be read by other programs.

    Args:
        chk_file_path (str | File): Path to the Gaussian checkpoint file or File object.

    Returns:
        str: Path to the generated formatted checkpoint file (.fchk).

    Raises:
        RuntimeError: If the formchk command fails.
        FileNotFoundError: If the checkpoint file doesn't exist.

    Note:
        The output .fchk file is created in the same directory as the input .chk file
        with the same base name but .fchk extension.
    """
    chk_file = File(chk_file_path)
    fchk_file_path = os.path.join(chk_file.file_dir, chk_file.file_prefix + ".fchk")
    shell_run(f"formchk {chk_file_path} > /dev/null")
    logger.debug(f"fchk file {fchk_file_path} is saved.")
    return fchk_file_path


def _get_atom_lines(mol2_block: str):
    """Extracts and parses atom lines from a MOL2 format block.

    Parses the atom section of a MOL2 format string and returns the atom
    information as a list of tokenized lines.

    Args:
        mol2_block (str): Complete MOL2 format string containing molecular data.

    Returns:
        List[List[str]]: List of atom lines where each line is split into tokens.
            Each atom line contains at least 9 fields including coordinates and charges.

    Raises:
        ValueError: If the MOL2 block doesn't contain required @<TRIPOS>ATOM or
            @<TRIPOS>BOND sections.

    Note:
        Only returns atom lines that have at least 9 fields to ensure complete
        atomic information including partial charges.
    """
    lines = mol2_block.splitlines()
    start_idx = lines.index("@<TRIPOS>ATOM")
    end_idx = lines.index("@<TRIPOS>BOND")
    atom_lines = lines[start_idx + 1 : end_idx]
    return [line.split() for line in atom_lines if len(line.split()) >= 9]


def _get_mol2_lines_with_charge(mol2_block: str, charge_list: list[float]):
    """Updates partial charges in a MOL2 format block with new charge values.

    Replaces the partial charges (9th column) in the atom section of a MOL2 format
    string with values from the provided charge list.

    Args:
        mol2_block (str): Complete MOL2 format string containing molecular data.
        charge_list (list[float]): List of new partial charges to assign to atoms. The list is
            consumed (modified) during processing via pop(0) operations.

    Returns:
        str: Modified MOL2 format string with updated partial charges.

    Raises:
        IndexError: If charge_list has fewer charges than atoms in the MOL2 block.
        ValueError: If the MOL2 block format is invalid.

    Note:
        The charge_list is modified during execution as charges are consumed.
        Only atoms with complete information (≥9 fields) have their charges updated.
        The formatting follows MOL2_LINE_FMT specification for proper alignment.
    """
    result_lines = []
    atom_line_flag = False
    for line in mol2_block.splitlines():
        if line.startswith("@<TRIPOS>ATOM"):
            atom_line_flag = True
            result_lines.append(line)
            continue
        elif line.startswith("@<TRIPOS>BOND"):
            atom_line_flag = False
        if atom_line_flag:
            line_contents = line.split()
            if len(line_contents) >= 9:
                line_contents[8] = str(charge_list.pop(0))
                line = MOL2_LINE_FMT.format(*line_contents)
        result_lines.append(line)
    return "\n".join(result_lines) + "\n"


def merge_mol2_charge(charge_mol2_block: str, origin_mol2_block: str) -> str:
    """Merges partial charges from one MOL2 block into another MOL2 block.

    Extracts partial charges from the charge_mol2_block and applies them to the
    corresponding atoms in the origin_mol2_block, preserving the original structure
    and connectivity while updating the charges.

    Args:
        charge_mol2_block (str): MOL2 format string containing the source charges.
        origin_mol2_block (str): MOL2 format string to receive the new charges.

    Returns:
        str: Modified MOL2 format string with charges from charge_mol2_block
            applied to the structure from origin_mol2_block.

    Raises:
        IndexError: If the number of atoms differs between the two MOL2 blocks.
        ValueError: If either MOL2 block has invalid format.

    Note:
        Both MOL2 blocks must have the same number of atoms with complete
        atomic information for proper charge transfer.
    """
    chg_list = [line[8] for line in _get_atom_lines(charge_mol2_block)]
    return _get_mol2_lines_with_charge(origin_mol2_block, chg_list)
