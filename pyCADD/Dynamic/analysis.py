"""Molecular dynamics trajectory analysis utilities.

This module provides functions for analyzing molecular dynamics trajectories including
RMSD/RMSF calculations, hydrogen bond analysis, distance and angle tracking, energy
analysis, and frame extraction capabilities.
"""

import logging
import os
from typing import Any

import pandas as pd
import pytraj as pt
from numpy import array, ndarray

from pyCADD.Dynamic.constant import PROCESS_MDOUT_PERL
from pyCADD.utils.common import BaseFile, ChDir
from pyCADD.utils.tool import shell_run, write_file

logger = logging.getLogger(__name__)


def _calc_rmsd(
    trajectory: Any, mask: str = "@CA", reference: Any = 0, save_dir: str = None, *args, **kwargs
) -> ndarray:
    """Calculates Root Mean Square Deviation (RMSD) for a molecular trajectory.

    Computes RMSD values for each frame of the trajectory after superposition.
    The trajectory is first aligned to a reference structure, then RMSD is calculated
    for the specified atoms.

    Args:
        trajectory (Any): PyTraj trajectory object containing molecular dynamics data.
        mask (str, optional): Atom selection mask for RMSD calculation. Defaults to "@CA" (alpha carbons).
        reference (Any, optional): Reference frame or structure for RMSD calculation. Defaults to 0 (first frame).
        save_dir (str, optional): Directory to save results. Defaults to current working directory.
        *args: Additional positional arguments passed to pytraj.rmsd.
        **kwargs: Additional keyword arguments passed to pytraj.rmsd.

    Returns:
        ndarray: Array of RMSD values for each frame in the trajectory.

    Note:
        Results are automatically saved as "RMSD_RESULTS.csv" in the specified directory.
        The trajectory is superposed using the same mask before RMSD calculation.
    """
    trajectory = pt.superpose(trajectory, mask=mask, ref=0)
    data_rmsd = pt.rmsd(trajectory, mask=mask, ref=reference, *args, **kwargs)
    save_dir = save_dir or os.getcwd()
    os.makedirs(save_dir, exist_ok=True)
    output_file = os.path.join(save_dir, "RMSD_RESULTS.csv")
    pd.DataFrame(data_rmsd, columns=["rmsd"]).to_csv(output_file, index=False)
    return data_rmsd


def _calc_rmsf(
    trajectory: Any,
    mask: str = "@CA",
    options: str = "byres",
    save_dir: str = None,
    *args,
    **kwargs,
) -> ndarray:
    """Calculates Root Mean Square Fluctuation (RMSF) for molecular trajectory.

    Computes RMSF values to measure the flexibility of atoms or residues over
    the course of the simulation. The trajectory is first aligned to remove
    global translation and rotation.

    Args:
        trajectory (Any): PyTraj trajectory object containing molecular dynamics data.
        mask (str, optional): Atom selection mask for RMSF calculation. Defaults to "@CA" (alpha carbons).
        options (str, optional): RMSF calculation options. Defaults to "byres" (per residue).
        save_dir (str, optional): Directory to save results. Defaults to current working directory.
        *args: Additional positional arguments passed to pytraj.rmsf.
        **kwargs: Additional keyword arguments passed to pytraj.rmsf.

    Returns:
        ndarray: Array of RMSF values, typically containing residue numbers and RMSF values.

    Note:
        Results are automatically saved as "RMSF_RESULTS.csv" with columns ["resnum", "rmsf"].
        Higher RMSF values indicate greater flexibility of the corresponding atoms/residues.
    """
    trajectory = pt.superpose(trajectory, mask=mask, ref=0)
    data_rmsf = pt.rmsf(trajectory, mask=mask, options=options, *args, **kwargs)
    save_dir = save_dir or os.getcwd()
    os.makedirs(save_dir, exist_ok=True)
    output_file = os.path.join(save_dir, "RMSF_RESULTS.csv")
    pd.DataFrame(data_rmsf, columns=["resnum", "rmsf"]).to_csv(output_file, index=False)
    return data_rmsf


def _calc_hbond(
    trajectory: Any,
    mask: str = ":*",
    distance: float = 3.0,
    angle: float = 135,
    options: str = None,
    save_dir: str = None,
    *args,
    **kwargs,
) -> tuple:
    """Analyzes hydrogen bonds in molecular dynamics trajectory.

    Performs comprehensive hydrogen bond analysis including identification of H-bonds,
    calculation of distances and angles, and generation of lifetime statistics using
    both PyTraj and cpptraj tools.

    Args:
        trajectory (Any): PyTraj trajectory object containing molecular dynamics data.
        mask (str, optional): Amber-style atom selection mask. Defaults to ":*" (all residues).
        distance (float, optional): Distance cutoff for H-bond detection in Angstroms. Defaults to 3.0.
        angle (float, optional): Angle cutoff for H-bond detection in degrees. Defaults to 135.
        options (str, optional): Additional H-bond analysis options. If None, uses default options.
        save_dir (str, optional): Directory to save results. Defaults to current working directory.
        *args: Additional positional arguments passed to pytraj.hbond.
        **kwargs: Additional keyword arguments passed to pytraj.hbond.

    Returns:
        tuple: A tuple containing (hbond_distance, hbond_angle) arrays with H-bond
            distance and angle data for each frame.

    Note:
        Generates multiple output files:
        - HBOND_RESULTS.dat: Average H-bond statistics
        - HBOND_NUM_TIME.csv: Number of H-bonds vs time
        - HBOND_DIS_RESULTS.csv: H-bond distances over time
        - HBOND_ANG_RESULTS.csv: H-bond angles over time
        - lifetime.dat: H-bond lifetime statistics from cpptraj
    """
    save_dir = save_dir or os.getcwd()
    os.makedirs(save_dir, exist_ok=True)
    avgout_file = f"{save_dir}/HBOND_RESULTS.dat"
    options = f"avgout {avgout_file} printatomnum nointramol" if options is None else options

    logger.info(f"Amber mask: {mask}")
    logger.info(f"Distance cutoff: {distance} angstrom")
    logger.info(f"Angle cutoff: {angle}")
    logger.info(f"Options: {options}")
    hbond = pt.hbond(
        trajectory, mask=mask, distance=distance, angle=angle, options=options, *args, **kwargs
    )

    # distance_mask from get_amber_mask() is acceptor-donorH
    # but needed distance is calculated from acceptor-donor
    # distance_mask, angle_mask = hbond.get_amber_mask()
    _, angle_mask = hbond.get_amber_mask()
    distance_mask = [(_mask.split()[0], _mask.split()[2]) for _mask in angle_mask]
    distance_mask = [" ".join(_mask_tuple) for _mask_tuple in distance_mask]
    distance_mask = array(distance_mask)

    hbond_distance = _trace_distance(trajectory, distance_mask)
    hbond_angle = _trace_angle(trajectory, angle_mask)

    hbond_num_time_file = os.path.join(save_dir, "HBOND_NUM_TIME.csv")
    hbond_distance_file = os.path.join(save_dir, "HBOND_DIS_RESULTS.csv")
    hbond_angle_file = os.path.join(save_dir, "HBOND_ANG_RESULTS.csv")

    # generate Hbond lifetime with cpptraj
    logger.info("Generating Hbond lifetime with cpptraj...")
    cpptraj_input_file = os.path.join(save_dir, "cpptraj.in")
    lifetime_file = os.path.join(save_dir, "lifetime.dat")
    cpptraj_input = f"parm {trajectory.topology.filename}\n"
    cpptraj_input += f"trajin {trajectory.filename}\n"
    cpptraj_input += f"hbond hbonds {mask} dist {distance} angle {angle} printatomnum nointramol series uuseries {lifetime_file}\n"
    cpptraj_input += f"run\nquit\n"
    write_file(cpptraj_input_file, cpptraj_input)
    shell_run(f"cpptraj -i {cpptraj_input_file}")

    pd.DataFrame(hbond.total_solute_hbonds()).to_csv(hbond_num_time_file)
    pd.DataFrame(hbond_distance, index=distance_mask).T.to_csv(hbond_distance_file, index=False)
    pd.DataFrame(hbond_angle, index=angle_mask).T.to_csv(hbond_angle_file, index=False)
    return hbond_distance, hbond_angle


def _trace_distance(
    trajectory: Any, mask: str, save: bool = False, save_dir: str = None, *args, **kwargs
) -> ndarray:
    """Traces distances between specified atoms throughout the trajectory.

    Calculates distances between pairs of atoms or groups of atoms for each frame
    in the molecular dynamics trajectory. Useful for monitoring specific
    interactions or conformational changes.

    Args:
        trajectory (Any): PyTraj trajectory object containing molecular dynamics data.
        mask (str): Amber-style atom selection mask defining the distance measurement.
            Format: "atom1 atom2" for distance between two atoms.
        save (bool, optional): Whether to save results to CSV file. Defaults to False.
        save_dir (str, optional): Directory to save results. Defaults to current working directory.
        *args: Additional positional arguments passed to pytraj.distance.
        **kwargs: Additional keyword arguments passed to pytraj.distance.

    Returns:
        ndarray: Array of distance values for each frame in the trajectory.

    Note:
        When save=True, results are saved as "TRACE_DIS_RESULTS_{mask}.csv".
        The mask format should follow Amber atom selection syntax.
    """
    save_dir = save_dir or os.getcwd()
    os.makedirs(save_dir, exist_ok=True)
    output_file = os.path.join(save_dir, f"TRACE_DIS_RESULTS_{mask}.csv")
    data_trace_distance = pt.distance(trajectory, mask=mask, *args, **kwargs)
    if save:
        pd.DataFrame(data_trace_distance, index=mask).T.to_csv(output_file, index=False)
    return data_trace_distance


def _trace_angle(
    trajectory: Any, mask: str, save: bool = False, save_dir: str = None, *args, **kwargs
) -> ndarray:
    """Traces angles between specified atoms throughout the trajectory.

    Calculates angles defined by three atoms for each frame in the molecular
    dynamics trajectory. Useful for monitoring conformational changes,
    bond angles, or dihedral angles.

    Args:
        trajectory (Any): PyTraj trajectory object containing molecular dynamics data.
        mask (str): Amber-style atom selection mask defining the angle measurement.
            Format: "atom1 atom2 atom3" for angle between three atoms.
        save (bool, optional): Whether to save results to CSV file. Defaults to False.
        save_dir (str, optional): Directory to save results. Defaults to current working directory.
        *args: Additional positional arguments passed to pytraj.angle.
        **kwargs: Additional keyword arguments passed to pytraj.angle.

    Returns:
        ndarray: Array of angle values (in degrees) for each frame in the trajectory.

    Note:
        When save=True, results are saved as "TRACE_ANG_RESULTS_{mask}.csv".
        The mask should specify three atoms that define the angle measurement.
    """
    save_dir = save_dir or os.getcwd()
    os.makedirs(save_dir, exist_ok=True)
    output_file = os.path.join(save_dir, f"TRACE_ANG_RESULTS_{mask}.csv")
    data_trace_angle = pt.angle(trajectory, mask=mask, *args, **kwargs)
    if save:
        pd.DataFrame(data_trace_angle, index=mask).T.to_csv(output_file, index=False)
    return data_trace_angle


def _get_lowest_energy_info(mdout_file: BaseFile, save_dir: str = None) -> tuple:
    """Identifies the frame with the lowest potential energy from MD output.

    Processes an Amber mdout file to extract energy information and identifies
    the simulation frame with the lowest total potential energy. This is useful
    for finding the most stable configuration during the simulation.

    Args:
        mdout_file (BaseFile): Amber mdout file containing energy data from MD simulation.
        save_dir (str, optional): Directory to save results. Defaults to current working directory.

    Returns:
        tuple: A tuple containing (frame_index, time, energy) where:
            - frame_index (int): Frame number with lowest energy
            - time (float): Simulation time at the lowest energy frame
            - energy (float): Lowest potential energy value

    Note:
        Uses a Perl script to process the mdout file and extract energy data.
        Results are saved as "LOWEST_ENERGY_ST_RESULTS.csv".
        If the first frame has the lowest energy, the function looks for the
        second-lowest to avoid initialization artifacts.

    Raises:
        FileNotFoundError: If the mdout file doesn't exist.
        ValueError: If energy data cannot be parsed from the mdout file.
    """
    save_dir = save_dir or os.getcwd()
    os.makedirs(save_dir, exist_ok=True)
    with ChDir(save_dir):
        process_script_path = os.path.join(save_dir, "process_mdout.perl")
        write_file(process_script_path, PROCESS_MDOUT_PERL)
        os.system(f"chmod +x {process_script_path} && {process_script_path} {mdout_file.file_path}")
    energy_file = os.path.join(save_dir, "summary.EPTOT")
    energy_df = pd.read_csv(energy_file, sep="\s+", header=None, names=["time", "energy"])
    frame = energy_df["energy"].idxmin()
    if frame == 1:
        frame = energy_df.loc[2:, "energy"].idxmin()
    time, energy = energy_df.iloc[frame, :].values

    output_file = os.path.join(save_dir, "LOWEST_ENERGY_ST_RESULTS.csv")
    result_str = f"frame,time,energy\n"
    result_str += f"{frame},{time},{energy}\n"
    write_file(output_file, result_str)
    return frame, time, energy


def _extract_frame(traj: pt.Trajectory, frame_indices: list[int], save_dir: str = None) -> list:
    """Extracts specific frames from a trajectory and saves them as PDB files.

    Extracts individual frames from a molecular dynamics trajectory and saves
    each frame as a separate PDB file. This is useful for structural analysis
    of specific conformations or time points.

    Args:
        traj (pt.Trajectory): PyTraj trajectory object to extract frames from.
        frame_indices (list[int]): List of frame indices to extract from the trajectory.
        save_dir (str, optional): Directory to save extracted PDB files. 
            Defaults to current working directory.

    Returns:
        list: List of file paths to the extracted PDB files.

    Note:
        Output files are named "FRAME_{frame_number}.pdb" where {frame_number}
        corresponds to the frame index in the trajectory. Existing files with
        the same names will be overwritten.

    Example:
        >>> frames = _extract_frame(traj, [0, 100, 200], "/output/dir")
        >>> # Creates: FRAME_0.pdb, FRAME_100.pdb, FRAME_200.pdb
    """
    save_dir = save_dir or os.getcwd()
    os.makedirs(save_dir, exist_ok=True)
    result_files = []
    output_file_template = os.path.join(save_dir, "FRAME_{frame}.pdb")
    for frame in frame_indices:
        current_output_file = output_file_template.format(frame=frame)
        pt.write_traj(
            current_output_file,
            traj=traj,
            frame_indices=[frame],
            overwrite=True,
        )
        result_files.append(current_output_file)
    return result_files
