import os
import signal
from typing import Tuple, Union

import pandas as pd

from pyCADD.utils.common import File, TimeoutError

from . import (
    Job,
    Structure,
    StructureReader,
    StructureWriter,
    center_of_mass,
    fileutils,
    jobcontrol,
    logger,
)
from .common import MaestroFile


def launch(cmd: str, timeout: int = 0) -> Job:
    """Launch a Schrodinger job.

    Args:
        cmd (str): command to be executed
        timeout (int, optional): timeout in seconds. Defaults to 0.

    Returns:
        schrodinger.job.Job: Schrodinger Job object
    """
    signal.signal(signal.SIGALRM, lambda signum, frame: TimeoutError)
    signal.alarm(timeout)
    cmd_list = cmd.split(" ")
    job = jobcontrol.launch_job(cmd_list)
    logger.debug(f"Command: {cmd}")
    logger.debug(f"JobId: {job.JobId}")
    logger.debug(f"Job Name: {job.Name}")
    logger.debug(f"Job Status: {job.Status}")
    try:
        job.wait()
    except TimeoutError:
        raise TimeoutError(f"Job {job.Name} timeout.")
    finally:
        signal.alarm(0)

    if not job.StructureOutputFile:
        logger.debug(f"Job failed: {job.Name}")
    else:
        logger.debug(f"File saved: {job.StructureOutputFile}")
    return job


def collect_structures(
    list_file: str, output_file: str = None, data_dir: str = None, required_file_ext: str = "maegz"
) -> None:
    """Collect and extract structures from data directory according to a list file including required structure names, and save them to another file.
    The structure names(without the extension) should be in the first column with a header.

    Args:
        list_file (str): list file path
        output_file (str, optional): output file path. Defaults to a .mae file with the same name as the list file.
        data_dir (str, optional): directory containing structure files. Defaults to current working directory.
        required_file_ext (str, optional): file extension of the required structure files. Defaults to 'maegz'.
    """
    data_dir = os.getcwd() if data_dir is None else data_dir
    list_file = File(list_file)
    output_file = File(output_file, exist=False)

    required_list = pd.read_csv(list_file.file_path).iloc[:, 0]
    required_sts = []

    for index, structure in enumerate(required_list):
        structure_name, structure_ext = os.path.splitext(os.path.basename(structure))
        required_file_ext = structure_ext.replace(".", "") if structure_ext else required_file_ext
        required_file_path = os.path.join(data_dir, f"{structure_name}.{required_file_ext}")
        if not os.path.exists(required_file_path):
            logger.warning(f"Structure file {required_file_path} not found!")
            continue
        st = StructureReader.read(required_file_path)
        st.property["i_user_SortedIndex"] = index
        required_sts.append(st)

    with StructureWriter(output_file.file_path) as writer:
        writer.extend(required_sts)


def convert_format(
    file_path: str,
    output_file: str,
    to_format: str = None,
    save_dir: str = None,
    overwrite: bool = False,
) -> str:
    """Convert the file to another format using Schrodinger API.

    Args:
        input_file_path (str): input file path
        output_file_path (str): output file path
        to_format (str): target format. If None, the format will be inferred from the output file extension.
        save_dir (str, optional): directory to save the converted file. Defaults to current working directory.
        overwrite (bool, optional): overwrite the existing file. Defaults to False.

    Raises:
        ValueError: Unsupported format
        FileNotFoundError: File to convert is not found
        FileExistsError: File already exists

    Returns:
        str: converted file path
    """
    _input_file = File(file_path)
    _output_file = File(output_file, exist=False)
    save_dir = _output_file.file_dir if save_dir is None else os.path.abspath(save_dir)
    to_format = _output_file.file_ext if to_format is None else to_format
    to_format = fileutils.get_structure_file_format(to_format)

    converted_file = os.path.join(save_dir, _output_file.file_name)
    if not overwrite and os.path.exists(converted_file):
        raise FileExistsError(f"File already exists: {converted_file}")
    try:
        st = StructureReader.read(_input_file.file_path)
        st.write(converted_file, format=to_format)
    except Exception:
        raise ValueError(f"Unsupported format: {to_format}")

    return converted_file


def get_centroid(
    file: Union[MaestroFile, Structure, str], structure_index: int = 0
) -> Tuple[float, float, float]:
    """Get the centroid of a structure file.

    Args:
        file_path (MaestroFile | str | Structure): file, maestro structure, or file path to read
        structure_index (int): Index of the structure to calculate the centroid. Defaults to 0.

    Returns:
        tuple[float, float, float]: centroid coordinates
    """
    if isinstance(file, str):
        file = MaestroFile(file)
    elif isinstance(file, Structure):
        return tuple(center_of_mass(file))
    return tuple(center_of_mass(file.structures[structure_index]))
