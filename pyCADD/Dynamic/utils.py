"""Utility functions for molecular dynamics simulation preparation and monitoring.

This module provides helper functions for charge calculation, parameter generation,
simulation progress tracking, and system analysis utilities used throughout
the molecular dynamics workflow.
"""

import logging
import os
import re
from time import sleep

from pyCADD.utils.common import File, ChDir
from pyCADD.utils.tool import _get_progress, shell_run

logger = logging.getLogger(__name__)


def calc_am1bcc(
    structure_file: str | File,
    charge: int,
    multiplicity: int = 1,
    force_field: str = "gaff2",
    save_dir: str = None,
) -> File:
    """Calculates AM1-BCC partial charges for a small molecule.

    Uses antechamber to calculate AM1-BCC charges and assigns GAFF2 atom types
    to a small molecule structure. This is a faster alternative to quantum
    mechanical charge calculation methods.

    Args:
        structure_file (str | File): Input molecular structure file (SDF, MOL2, PDB, etc.).
        charge (int): Net charge of the molecule.
        multiplicity (int, optional): Spin multiplicity of the molecule. Defaults to 1.
        force_field (str, optional): Force field to assign atom types. Defaults to "gaff2".
        save_dir (str, optional): Directory to save the output MOL2 file.
            Defaults to current working directory.

    Returns:
        File: MOL2 file with AM1-BCC charges and GAFF2 atom types.

    Raises:
        FileNotFoundError: If antechamber calculation fails or output file not generated.
        RuntimeError: If antechamber execution fails.

    Note:
        AM1-BCC charges are semi-empirical and faster than ab initio methods
        but may be less accurate than RESP charges for some applications.
        The output file is named "{input_prefix}_bcc.mol2".
    """
    save_dir = save_dir or os.getcwd()
    structure_file = File(structure_file)
    file_path = structure_file.file_path
    file_prefix = structure_file.file_prefix
    file_suffix = structure_file.file_suffix
    mol2_file_path = os.path.join(save_dir, f"{file_prefix}_bcc.mol2")
    logger.info(f"Starting AM1/BCC calculation for {structure_file.file_path}")
    with ChDir(save_dir):
        shell_run(
            f"antechamber -fi {file_suffix} -i {file_path} -fo mol2 -o {mol2_file_path} "
            f"-c bcc -nc {charge} -m {multiplicity} -at {force_field}"
        )
    if not os.path.exists(mol2_file_path):
        raise FileNotFoundError(f"AM1/BCC calculation failed, {mol2_file_path} not found.")
    return File(mol2_file_path)


def _get_water_resnum(pdb_file: File | str) -> list:
    """Extracts residue numbers of water molecules from solvated PDB file.

    Parses a solvated system PDB file to identify and return the residue numbers
    of all water molecules (WAT residues). This is useful for defining water
    masks in MD analysis.

    Args:
        pdb_file (File | str): Solvated system PDB file containing water molecules.

    Returns:
        list: List of residue numbers (as strings) for all water molecules in the system.

    Note:
        Uses grep and awk commands to extract residue numbers from the 5th column
        of PDB records containing "WAT" residue names. Assumes standard PDB format.

    Example:
        >>> water_nums = _get_water_resnum("system_solvated.pdb")
        >>> print(f"Found {len(water_nums)} water molecules")
    """
    pdb_file = File(pdb_file)
    return (
        os.popen("cat %s | grep WAT -w | awk '{print $5}'" % pdb_file.file_path).read().splitlines()
    )


def run_parmchk2(
    mol2_file_path: str | File,
    save_dir: str = None,
) -> File:
    """Generates missing force field parameters using parmchk2.

    Runs Amber's parmchk2 program to identify and generate missing force field
    parameters for a molecule based on its MOL2 file with assigned atom types.
    Creates a frcmod file containing additional parameters not found in the
    standard GAFF force field.

    Args:
        mol2_file_path (str | File): File or Path to MOL2 file with assigned atom types (from antechamber).
        save_dir (str, optional): Directory to save the frcmod parameter file.
            Defaults to current working directory.

    Returns:
        File: Generated frcmod file containing missing force field parameters.

    Raises:
        RuntimeError: If parmchk2 execution fails.
        FileNotFoundError: If input MOL2 file doesn't exist.

    Note:
        The frcmod file contains additional parameters (bonds, angles, dihedrals)
        that are not available in the standard GAFF parameter set. This file
        is used together with GAFF in tLEaP to build complete molecular systems.
        Output file is named "{mol2_prefix}.frcmod".
    """
    mol2_file = File(mol2_file_path)
    save_dir = save_dir or os.getcwd()
    file_prefix = mol2_file.file_prefix
    frcmod_file_path = os.path.join(save_dir, file_prefix + ".frcmod")
    cmd = f"parmchk2 -i {mol2_file.file_path} -f mol2 -o {frcmod_file_path}"
    shell_run(cmd)
    return File(frcmod_file_path)


def _trace_progress(output_file_path: str, step: int):
    """Tracks and displays progress of a running MD simulation.

    Monitors a molecular dynamics simulation by parsing the output file to extract
    the current step number and displays a real-time progress bar. This function
    runs continuously until the simulation completes.

    Args:
        output_file_path (str): Path to the MD simulation output file (.out).
        step (int): Total number of simulation steps expected.

    Note:
        The function monitors the output file by:
        - Looking for "NSTEP" lines to get current step number
        - Looking for "Final" keyword to detect completion
        - Updating a progress bar every second
        - Automatically stopping when simulation finishes

        This is typically used with background (nohup) simulations to provide
        real-time feedback on simulation progress.

    Raises:
        FileNotFoundError: If the output file doesn't exist or isn't readable.
    """
    progress, taskID = _get_progress("Molecule Dynamics Simulation", "bold cyan", total=step)
    progress.start()
    progress.start_task(taskID)
    sleep(5)
    while not progress.finished:
        _current = os.popen(f"tail -n 10 {output_file_path} | grep NSTEP").read()
        _finished = os.popen(f"tail -n 20 {output_file_path} | grep Final").read()
        if _current:
            current_step = re.findall(r"\d+", _current)[0]
            progress.update(taskID, completed=int(current_step))
        elif _finished:
            progress.update(taskID, completed=step)
            break
        sleep(1)
    progress.stop()
