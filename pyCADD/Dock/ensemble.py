import logging
import os
import multiprocessing
from shutil import move

# For Schrodinger 2021-2 or newer release
import importlib
from tkinter import N
from typing import Iterable
importlib.reload(multiprocessing)

from pyCADD.utils.tool import download_pdb, download_pdb_list, makedirs_from_list, _get_progress, get_config
from pyCADD.Dock.common import PDBFile, MaestroFile, GridFile, DockResultFile, ComplexFile, ReceptorFile, LigandFile, MultiInputFile
from pyCADD.Dock.core import minimize, grid_generate, dock, calc_admet, calc_mmgbsa, calc_volume
from pyCADD.Dock.data import extra_docking_data

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
    
def creat_mapping(grid_file_list:list, ligand_file_list:list) -> tuple:
    '''
    将所有受体与全部配体小分子建立完全映射关系

    Parameters
    ----------
    grid_file_list : list
        受体格点文件列表
    ligand_file_list : list
        配体文件列表
    '''
    mapping_results = []
    logger.debug(f'Prepare to map {len(ligand_file_list)} ligands to {len(grid_file_list)} receptors')

    for grid_file in grid_file_list:
        mapping_results.extend([(grid_file, ligand_file) for ligand_file in ligand_file_list])

    return mapping_results

def _multiprocssing_run(func, _iterable:Iterable, *args, job_name:str, num_parallel:int):
    '''
    多进程运行函数 并自动生成进度条

    Parameters
    ----------
    func : function
        运行函数
    _iterable : Iterable
        可迭代对象 即函数作用对象总集
    *args
        传入func函数的其他参数
    job_name : str
        进程名称
    num_parallel : int
        进程数量
    
    Returns
    -------
    List
        所有成功完成任务的返回值
    
    '''
    progress, taskID = _get_progress(job_name, 'bold cyan', len(_iterable))
    returns = []
    def success_handler(result):
        '''
        处理成功的任务
        '''
        returns.append(result)
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
    return returns
 
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

        self.pdbfile_list = None
        self.minimized_file_list = None
        self.ligand_list = None
        self.mapping = None
        self.grid_file_list = None
        self.ligand_file_list = None
        self.dock_file_list = None

        self.pdb_save_dir = os.path.join(os.getcwd(), 'pdb')
        self.minimize_save_dir = os.path.join(os.getcwd(), 'minimize')
        self.grid_save_dir = os.path.join(os.getcwd(), 'grid')
        self.ligand_save_dir = os.path.join(os.getcwd(), 'ligands')
        self.complex_save_dir = os.path.join(os.getcwd(), 'complex')
        self.protein_save_dir = os.path.join(os.getcwd(), 'protein')
        self.base_dock_save_dir = os.path.join(os.getcwd(), 'dockfiles')
        self.result_save_dir = os.path.join(os.getcwd(), 'result')
        
        makedirs_from_list([
            self.pdb_save_dir, 
            self.minimize_save_dir, 
            self.grid_save_dir,
            self.ligand_save_dir,
            self.complex_save_dir,
            self.protein_save_dir,
            self.base_dock_save_dir,
            self.result_save_dir
            ])
    
    def download_all_pdb(self, overwrite:bool=False) -> None:
        '''
        下载列表中的所有PDB文件

        Parameters
        ----------
        overwrite : bool
            是否覆盖已存在的文件
        '''
        pdbid_list = self.pdbid_list
        pdb_save_dir = self.pdb_save_dir
        download_pdb_list(pdbid_list, pdb_save_dir, overwrite)
        self.pdbfile_list = [PDBFile(pdbfile_path) for pdbfile_path in self.input_file.get_pdbfile_path_list(pdb_save_dir)]

    def keep_single_chain(self):
        '''
        将所有PDB文件转换为单链构象
        '''
        logger.debug(f'Prepare to keep single chain for {len(self.pdbid_list)} PDB files')
        self.pdbfile_list = [PDBFile(pdbfile_path).keep_single_chain() for pdbfile_path in self.input_file.get_pdbfile_path_list(self.pdb_save_dir)]

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
        pairs_list = self.pairs_list
        minimize_save_dir = self.minimize_save_dir

        self.download_all_pdb(overwrite)
        pdbfile_list = self.pdbfile_list
        if keep_single_chain:
            self.keep_single_chain()

        self.minimized_file_list = _multiprocssing_run(minimize, pdbfile_list, side_chain, missing_loop, del_water, minimize_save_dir, overwrite, job_name='Minimizing Structures', num_parallel=num_parallel)
        # self.minimized_list = [os.path.join(minimize_save_dir, f'{pdbid}_minimized.mae') for pdbid in pdbid_list]
        # 某一结构优化失败时将抛出FileNotFoundError
        # self.minimized_file_list = [ComplexFile(os.path.join(minimize_save_dir, f'{pdbid}_minimized.mae'), ligand) for pdbid, ligand in pairs_list]
        return self.minimized_file_list

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
            是否覆盖已存在的zip文件
        
        Returns
        -------
        list
            生成的网格文件列表
        '''
        pairs_list = self.pairs_list
        logger.debug(f'Prepare to generate grids for {len(pairs_list)} structures')

        _pairs_list = self.minimized_file_list
        grid_save_dir = self.grid_save_dir

        if _pairs_list is None:
            raise FileNotFoundError('No minimized file found. Please run minimize first.')

        self.grid_file_list = _multiprocssing_run(grid_generate, _pairs_list, gridbox_size, grid_save_dir, overwrite, job_name='Generating Grids', num_parallel=num_parallel)
        # self.grid_list = self.input_file.get_gridfile_path_list(grid_save_dir)
        
        return self.grid_file_list

    def minimized_split(self) -> None:
        '''
        将优化的结构拆分为ligand和protein 保存至相应位置
        '''
        minimized_file_list = self.minimized_file_list
        if minimized_file_list is None:
            raise RuntimeError('No minimized file found. Please run minimize first.')

        logger.debug(f'Prepare to split {len(minimized_file_list)} structures')
        logger.debug(f'Protein file will be saved in {self.protein_save_dir}')
        logger.debug(f'Ligand file will be saved in {self.ligand_save_dir}')
        logger.debug(f'Complex file will be saved in {self.complex_save_dir}')
        # split不进行多进程化
        for minimized_file in minimized_file_list:
            minimized_file.split(protein_dir=self.protein_save_dir, ligand_dir=self.ligand_save_dir, complex_dir=self.complex_save_dir)
    
    def ligand_split(self, external_ligand_file:str, overwrite:bool=False) -> list:
        '''
        拆分外部ligand文件
        
        Parameter
        ---------
        external_ligand_file : str
            外部ligand文件路径
        
        Returns
        -------
        list
            拆分后的ligand文件列表
        '''
        logger.debug(f'Prepare to split {external_ligand_file}')
        ligand_save_dir = self.ligand_save_dir
        self.ligand_list = split_ligand(external_ligand_file, ligand_save_dir, overwrite)
        return self.ligand_list
        
    def creat_mapping(self):
        '''
        建立映射
        '''
        if self.ligand_list is None:
            raise RuntimeError('Please run ligand_split first.')
        
        self.grid_file_list = [GridFile(os.path.join(self.grid_save_dir, f'{pdbid}_glide-grid_{ligname}.zip') for pdbid, ligname in self.pairs_list)]
        self.ligand_file_list = [LigandFile(os.path.join(self.ligand_save_dir, f'{ligand}.mae')) for ligand in self.ligand_list]
        self.mapping = creat_mapping(self.grid_file_list, self.ligand_file_list)

    def multi_dock(self, precision:str='SP', calc_rmsd:bool=False, num_parallel:int=NUM_PARALLEL, overwrite:bool=False) -> list:
        '''
        使用多进程调用Glide 执行批量分子对接

        Parameter
        ---------
        precision : str
            分子对接精度
        calc_rmsd : bool
            是否计算rmsd
        num_parallel : int
            并行进程数
        overwrite : bool
            是否覆盖已存在的mae文件
        '''
        logger.debug(f'Grid files in {self.grid_save_dir} will be used.')
        logger.debug(f'Ligand files in {self.ligand_save_dir} will be used.')
        logger.debug(f'Number of all jobs: {len(self.grid_file_list) * len(self.ligand_file_list)}')
        logger.debug(f'Docking precision: {precision}')
        logger.debug(f'Calculate rmsd: {calc_rmsd}')
        logger.debug(f'Number of parallel jobs: {num_parallel}')

        if self.mapping is None:
            raise RuntimeError('Please run creat_mapping first.')
        
        self.dock_file_list = _multiprocssing_run(dock, self.mapping, precision, calc_rmsd, self.base_dock_save_dir, overwrite, job_name='Ensemble Docking', num_parallel=num_parallel)
        
    def multi_extract_data(self, precision:str='SP', num_parallel:int=NUM_PARALLEL, overwrite:bool=False) -> list:
        '''
        多进程 提取对接结果数据
        '''
        if self.dock_file_list is None:
            self.dock_file_list = [DockResultFile(os.path.join(self.base_dock_save_dir, pdbid, f'{pdbid}_{ligid}_glide-dock_{ligand}_{precision}.maegz')) for pdbid, ligid in self.pairs_list for ligand in self.ligand_list]
        
        total_data_list = _multiprocssing_run(extra_docking_data, self.dock_file_list, job_name='Extract Docking Data', num_parallel=NUM_PARALLEL)
        
        return total_data_list
