import logging
import os
import multiprocessing
from shutil import move

# For Schrodinger 2021-2 or newer release
import importlib
from typing import Iterable
importlib.reload(multiprocessing)

from pyCADD.utils.tool import download_pdb, download_pdb_list, makedirs_from_list, _get_progress, get_config
from pyCADD.Dock.common import PDBFile, MaestroFile, GridFile, DockResultFile, ComplexFile, ReceptorFile, LigandFile, MultiInputFile
from pyCADD.Dock.core import keep_chain, minimize, grid_generate, dock, calc_admet, calc_mmgbsa, calc_volume

logger = logging.getLogger(__name__)
NUM_PARALLEL = multiprocessing.cpu_count() // 4 * 3

def split_ligand(ligand_file:LigandFile, save_dir:str=None, overwrite:bool=False) -> list:
    '''
    将单个maestro文件中包含的所有小分子拆分为多个独立mae文件
    
    Parameters
    ----------
    ligand_file : LigandFile
        单个maestro文件
    save_dir : str
        拆分后的mae文件保存的目录
    overwrite : bool
        是否覆盖已存在的mae文件
    
    Returns
    -------
    list
        拆分后的配体名称列表   
            配体名称由 唯一索引index + ligand_name 组成 
    '''

    save_dir = ligand_file.file_dir if save_dir is None else save_dir
    logger.debug(f'Prepare to split structure file {ligand_file.file_name} to {save_dir}')

    progress, taskID = _get_progress('Reading Ligands', 'bold cyan', len(ligand_file.structures))
    progress.start()
    progress.start_task(taskID)

    label_list = []
    ligand_list = []
    activity_label_name = [f'{_type}_user_{_label}' for _type in ('b', 's') for _label in ('Activity', 'activity')]
    for index, structure in enumerate(ligand_file.structures):
        st_name = f"{index}-{structure.property['s_m_title']}"
        structure.property['i_user_StructureIndex'] = index

        st_activity = ''
        for _label in activity_label_name:
            try:
                st_activity = structure.property[_label]
                break
            except KeyError:
                continue
        label_list.append(f'{st_name},{st_activity}')

        output_file = os.path.join(save_dir, f'{st_name}.mae')
        if not os.path.exists(output_file) or overwrite:
            structure.write(output_file)
        
        ligand_list.append(st_name)
        progress.update(taskID, advance=1)

    with open(os.path.join(save_dir, 'label.csv'), 'w') as f:
        f.write('\n'.join(label_list))
    
    progress.stop()

    return ligand_list
    
def creat_mapping(pairs_list:list, ligand_list:list) -> tuple:
    '''
    将所有受体与全部配体小分子建立完全映射关系

    Parameters
    ----------
    pairs_list : list
        受体信息列表
    ligand_list : list
        配体名称列表
    '''
    mapping_results = []
    logger.debug(f'Prepare to map {len(ligand_list)} ligands to {len(pairs_list)} receptors')

    for pdbid, ligid in pairs_list:
        mapping_results.extend([(pdbid, ligid, ligand) for ligand in ligand_list])

    return mapping_results

def _multiprocssing_run(job_name:str, func, _iterable:Iterable, num_parallel:int=NUM_PARALLEL, *args):
    '''
    多进程运行函数 并自动生成进度条

    Parameters
    ----------
    job_name : str
        进程名称
    func : function
        运行函数
    _iterable : Iterable
        可迭代对象 即函数作用对象总集
    num_parallel : int
        进程数量
    *args
        传入func函数的其他参数

    '''
    progress, taskID = _get_progress(job_name, 'bold cyan', len(_iterable))

    def success_handler(result):
        '''
        处理成功的任务
        '''
        progress.update(taskID, advance=1)

    def _error_handler(error:Exception):
        '''
        异常处理函数
        '''
        logger.error(f'{error}')
        progress.update(taskID, advance=1)

    progress.start()
    progress.start_task(taskID)

    pool = multiprocessing.Pool(num_parallel, maxtasksperchild=1)
    for item in _iterable:
        if isinstance(item, Iterable):
            pool.apply_async(func, args=(*item, *args), callback=success_handler, error_callback=_error_handler)
        else:
            pool.apply_async(func, args=(item, *args), callback=success_handler, error_callback=_error_handler)
    pool.close()
    pool.join()

    progress.stop()
 
class Console:
    '''
    Ensemble docking 控制台对象
    '''
    def __init__(self, input_file:MultiInputFile) -> None:
        self.input_file = input_file
        self.pairs_list = input_file.get_pairs_list()
        self.pdbid_list = input_file.get_pdbid_list()
        self.grid_list = None
        self.minimized_list = None

        self.pdb_save_dir = os.path.join(os.getcwd(), 'pdb')
        self.minimize_save_dir = os.path.join(os.getcwd(), 'minimize')
        self.grid_save_dir = os.path.join(os.getcwd(), 'grid')
        self.ligand_save_dir = os.path.join(os.getcwd(), 'ligands')
        self.complex_save_dir = os.path.join(os.getcwd(), 'complex')
        self.protein_save_dir = os.path.join(os.getcwd(), 'protein')
        
        makedirs_from_list([
            self.pdb_save_dir, 
            self.minimize_save_dir, 
            self.grid_save_dir,
            self.ligand_save_dir,
            self.complex_save_dir,
            self.protein_save_dir
            ])

    def multi_minimize(self, keep_single_chain:bool=True, num_parallel:int=NUM_PARALLEL, side_chain:bool=True, missing_loop:bool=True, del_water:bool=True, overwrite:bool=False) -> list:
        '''
        使用多进程调用prepwizard 运行多个受体结构的优化
        
        Parameter
        ---------
        keep_single_chain : bool
            是否保留单链
        num_parallel : int
            并行进程数
        overwrite : bool
            是否覆盖已存在的mae文件
        
        Returns
        -------
        list
            优化后的mae文件列表
        '''
        logger.debug(f'Prepare to optimize and minimize {len(self.pairs_list)} structures')

        pdbid_list = self.pdbid_list
        pdb_save_dir = self.pdb_save_dir
        minimize_save_dir = self.minimize_save_dir

        download_pdb_list(pdbid_list, pdb_save_dir, overwrite)

        pdbfile_list = [PDBFile(pdbfile_path) for pdbfile_path in self.input_file.get_pdbfile_path_list(pdb_save_dir)]
        if keep_single_chain:
            pdbfile_list = [pdbfile.keep_chain() for pdbfile in pdbfile_list]

        _multiprocssing_run('Minimizing Structures', minimize, pdbfile_list, num_parallel, side_chain, missing_loop, del_water, minimize_save_dir, overwrite)
        self.minimized_list = [os.path.join(minimize_save_dir, f'{pdbid}_minimized.mae') for pdbid in pdbid_list]
        
        return self.minimized_list

    def multi_grid_generate(self, gridbox_size:int=20, num_parallel:int=NUM_PARALLEL, overwrite:bool=False) -> list:
        '''
        使用多进程调用Glide 运行多个受体结构的格点文件生成
        
        Parameter
        ---------
        gridbox_size : int
            格点大小
        num_parallel : int
            并行进程数
        overwrite : bool
            是否覆盖已存在的mae文件
        
        Returns
        -------
        list
            生成的网格文件列表
        '''
        pairs_list = self.pairs_list
        logger.debug(f'Prepare to generate grids for {len(pairs_list)} structures')

        minimized_file_save_dir = self.minimize_save_dir
        grid_save_dir = self.grid_save_dir

        try:
            _pairs_list = [(ComplexFile(os.path.join(minimized_file_save_dir, f'{pdbid}_minimized.mae')), ligid) for pdbid, ligid in pairs_list]
        except FileNotFoundError:
            raise FileNotFoundError('No minimized file found. Please run minimize first.')

        _multiprocssing_run('Generating Grids', grid_generate, _pairs_list, num_parallel, gridbox_size, grid_save_dir, overwrite)
        self.grid_list = self.input_file.get_gridfile_path_list(grid_save_dir)
        
        return self.grid_list

    def multi_split(self):
        pairs_list = self.pairs_list
        if self.minimized_list is None:
            raise ValueError('No minimized file found. Please run minimize first.')
        minimized_list = self.minimized_list

        logger.debug(f'Prepare to split {len(pairs_list)} structures')
        logger.debug(f'Protein file will be saved in {self.protein_save_dir}')
        logger.debug(f'Ligand file will be saved in {self.ligand_save_dir}')

        # TODO