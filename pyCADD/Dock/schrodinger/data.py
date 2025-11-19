import os
from typing import List, Union

import pandas as pd

from pyCADD.Dock.schrodinger.common import DockResultFile
from pyCADD.Dock.schrodinger.config import DataConfig

from . import logger


def extract_docking_data(
    docking_result_file: Union[DockResultFile, str], data_config: DataConfig = None
) -> List[dict]:
    """Extract docking data from a docking result file.

    Args:
        docking_result_file (Union[DockResultFile, str]): docking result file or path to the file.
        data_config (DataConfig, optional): data extracting config. By default, the configuration corresponding to the docking precision is used.

    Returns:
        List[dict]: docking result data list.
    """
    if isinstance(docking_result_file, str):
        docking_result_file = DockResultFile(docking_result_file)
    if data_config is None:
        data_config = DataConfig(precision=docking_result_file.metadata.precision)

    logger.debug(f"Prepare to extract data from file {docking_result_file.file_path}")
    include_receptor = docking_result_file.metadata.get("include_receptor", False)
    raw_data = docking_result_file.get_raw_results()
    if raw_data:
        raw_data = raw_data[1:] if include_receptor else raw_data
    result_data_list = []
    data_dict = {
        "pdbid": docking_result_file.metadata.pdbid,
        "precision": docking_result_file.metadata.precision,
        "internal_ligand_name": docking_result_file.metadata.internal_ligand_name,
        "docking_ligand_name": docking_result_file.metadata.docking_ligand_name,
    }
    data = data_dict.copy()
    for raw_data_dict in raw_data:
        data.update({key: raw_data_dict.get(value, None) for key, value in data_config.items()})
    result_data_list.append(data)

    return result_data_list


def save_docking_data(
    docking_result_file: Union[DockResultFile, str],
    data_config: DataConfig = None,
    save_dir: str = None,
    overwrite: bool = False,
) -> None:
    """Save docking data to a CSV file. The file name is the same as the docking result file.

    Args:
        docking_result_file (Union[DockResultFile, str]): docking result file or path to the file.
        data_config (DataConfig, optional): data extracting config. By default, the configuration corresponding to the docking precision is used.
        save_dir (str, optional): directory to save the data. Defaults to None.
        overwrite (bool, optional): whether to overwrite the existing file. Defaults to False.

    Raises:
        FileExistsError: raise if the output file already exists.
    """
    if isinstance(docking_result_file, str):
        docking_result_file = DockResultFile(docking_result_file)
    data_list = extract_docking_data(docking_result_file, data_config)
    df = pd.DataFrame(data_list)
    save_dir = os.path.abspath(save_dir) if save_dir is not None else os.getcwd()
    output_file = os.path.join(save_dir, f"{docking_result_file.file_prefix}.csv")
    if os.path.exists(output_file) and not overwrite:
        raise FileExistsError(f"File already exists: {output_file}")
    df.to_csv(output_file, index=False)
