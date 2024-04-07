import os
import logging
import shutil
from typing import Union, Literal

from pyCADD.utils.common import File, ChDir
from pyCADD.utils.tool import shell_run
from pyCADD.Dock.schrodinger.common import MaestroFile, DockResultFile
from pyCADD.Dock.schrodinger.config import FORCE_FILED_DICT
from pyCADD.Dock.schrodinger.utils import launch

def minimize(
    file:Union[MaestroFile, str], 
    ph: float=7.4,
    force_field:Literal['OPLS4', 'OPLS3e', 'OPLS3', 'OPLS2005']='OPLS4',
    fill_side_chain:bool=True, 
    add_missing_loop:bool=True, 
    del_water:bool=True, 
    watdist:float=5.0,
    rmsd_cutoff:float=0.3,
    save_dir:str=None, 
    overwrite:bool=False
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
        file = File(path=file)
    save_dir = save_dir if save_dir else os.getcwd()
    force_field = FORCE_FILED_DICT[force_field]
    logging.debug(f'Prepare to minimize {file.file_path}')
    minimized_file = os.path.join(save_dir, f'{file.file_prefix}__minimized.mae')
    if not overwrite and os.path.exists(minimized_file):
        logging.debug(f'File already existed: {minimized_file}')
        return MaestroFile(minimized_file)
    
    with ChDir(save_dir):
        job_name = f'{file.file_prefix}-Minimize'
        prepwizard_command = f'prepwizard -f {force_field} -r {rmsd_cutoff} -propka_pH {ph} -disulfides -s -JOBNAME {job_name}'
        if fill_side_chain:
            prepwizard_command += ' -fillsidechains'
        if add_missing_loop:
            prepwizard_command += ' -fillloops'
        if del_water:
            prepwizard_command += f' -watdist {watdist}'
        else:
            prepwizard_command += '-keepfarwat'
        _minimized_file = f'{file.file_prefix}_minimized.mae'
        prepwizard_command += f' {file.file_path} {_minimized_file}'
        
        launch(prepwizard_command)
        
        try:
            shutil.move(_minimized_file, minimized_file)
        except FileNotFoundError:
            raise RuntimeError(f'Minimization Process Failed: {file.file_prefix}')
        else:
            logging.debug(f'Minimized file saved: {minimized_file}')
            return MaestroFile(minimized_file)