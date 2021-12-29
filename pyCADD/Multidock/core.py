import logging
import multiprocessing
import os

from pyCADD.Dock.core import cal_admet, dock, grid_generate, minimize, split_com, cal_mmgbsa
from pyCADD.Multidock.prepare import minimize_prepare, split_ligand
from pyCADD.utils.getinfo import get_pdbfile_path_list, get_project_dir
from pyCADD.utils.tool import mkdirs, _get_progress

from concurrent.futures import ProcessPoolExecutor

logger = logging.getLogger('pyCADD.Multidock.core')

def read_receptors(list_file_path:str) -> list:
    '''
    处理受体列表文件

    Parameter
    ---------
    list_file_path : str
        列表文件PATH 文件中含有多行的 逗号分隔的PDBID与配体ID
        example:
            3A9E,REA
            3KMZ,EQO
            3KMR,EQN
            4DQM,LUF
            1DKF,BMS
            5K13,6Q7
    
    Return
    ---------
    list
        从列表文件中读取的多组受体信息

    '''
    receptor_list = []
    logger.debug('Reading receptors file: %s' % list_file_path)

    with open(list_file_path, 'r') as f:
        # 读取PDB列表中的每一行PDBID与配体名称信息
        pdbs_withlig = f.read().splitlines()        

    for i in pdbs_withlig:                                      # 按逗号分割解析列表中的PDB ID与配体名称

        pdb = i.split(',')[0].strip().upper()                   # PDB ID
        lig = i.split(',')[1].strip().upper()                   # 配体名称
        receptor_list.append((pdb, lig))                        # 将每一对作为元组储存至列表receptor_list
    
    # 生成标签文件
    cwd = get_project_dir()
    with open(cwd + '/ligands/original.csv','w') as f:
        f.write('Ligand,activity\n')
        for pdbid, lig in receptor_list:
            f.write('%s-lig-%s,origin\n' % (pdbid, lig))

    return receptor_list

def read_ligands(maefile:str, dirname:str='./') -> list:
    '''
    处理配体 拆分mae文件为单个化合物的结构

    Parameter
    ---------
    maefile : str
        包含所有配体小分子的单个mae文件PATH
    dirname : str
        拆分的mae储存目录PATH

    Return
    ---------
    list
        多组小分子配体名称
    '''
    logger.debug('Loading ligands file: %s' % maefile)
    return split_ligand(maefile, dirname)


def map(receptor_list:list, ligand_list:list) -> tuple:
    '''
    将所有受体与全部配体小分子建立1:1独立映射关系
    Parameters
    ----------
    receptor_list : list
        受体信息列表
    ligand_list : list
        配体名称列表
    
    Return
    ----------
    tuple
        映射关系元组
    '''
    tup = ()
    for receptor, lig in receptor_list:
        ligand_list.append('%s-lig-%s' % (receptor, lig))

    logger.debug('There is/are %s receptor(s) to be mapped' % len(receptor_list))
    logger.debug('There is/are %s ligand(s) to be mapped' % len(ligand_list))
    logger.info('A map of %s X %s will be created.' % (len(receptor_list), len(ligand_list)))
    progress,task = _get_progress('Creating Map', 'bold cyan', len(receptor_list))
    progress.start()
    progress.start_task(task)
    def _update(*arg):
        progress.update(task, advance=1)

    for receptor, lig in receptor_list:
        for ligand in ligand_list:
            tup += ((receptor, lig, ligand),)
        _update()
    progress.stop()
    return tup
    

def multi_minimize(pdblist:list):
    '''
    多进程 处理多个受体结构的优化
    
    Parameter
    ---------
    pdblist : list
        受体PDB ID列表
    '''
    progress, task = _get_progress('Optimizing structure', 'bold cyan', len(pdblist))
    minimize_prepare(pdblist)
    progress.start()
    def _update(*arg):
        progress.update(task, advance=1)

    cwd = get_project_dir()
    logger.debug('Current working directory: %s' % cwd)
    # minimized文件存放目录
    minimize_dir = cwd + '/minimize/'
    logger.debug('Minimized files will be saved in %s' % minimize_dir)

    pdbfiles = get_pdbfile_path_list(pdblist)

    # 暂时进入minimize文件存放目录 以优化晶体并保存结构文件于此处
    os.chdir(minimize_dir)
    # 最大进程数为CPU核心数量 1:1
    pool = multiprocessing.Pool(os.cpu_count(), maxtasksperchild=1)
    # 进度条
    progress.start_task(task)
    for pdbfile in pdbfiles:
        pool.apply_async(minimize, (pdbfile,), error_callback=error_handler, callback=_update)

    pool.close()
    pool.join()
    progress.stop()
    # 返回原始工作目录
    os.chdir(cwd)
    

def multi_grid_generate(receptor_list:list):
    '''
    多进程 生成多个受体结构的格点文件

    Parameter
    ----------
    receptor_list : list
        受体信息列表(PDBID, 配体ID)
    '''

    cwd = get_project_dir()
    logger.debug('Current working directory: %s' % cwd)
    grid_dir = cwd + '/grid/'
    minimize_dir = cwd + '/minimize/'
    logger.debug('Grid files will be saved in %s' % grid_dir)
    logger.debug('Minimized files will be used in %s' % minimize_dir)

    progress, task = _get_progress('Grid generating', 'bold cyan', len(receptor_list))
    progress.start()
    progress.start_task(task)
    def _update(*arg):
        progress.update(task, advance=1)

    # 暂时进入grid文件存放目录 计算格点文件并保存结构于此处
    os.chdir(grid_dir)
    pool = multiprocessing.Pool(os.cpu_count(), maxtasksperchild=1)

    for pdbid, lig in receptor_list:
        st_file_path = minimize_dir + '%s_minimized.mae' % pdbid
        # 默认格点大小20A
        pool.apply_async(grid_generate, (pdbid, lig, st_file_path, 20), error_callback=error_handler, callback=_update)
    
    pool.close()
    pool.join()
    progress.stop()
    # 返回原始工作目录
    os.chdir(cwd)

def multi_split(receptor_list:list):
    '''
    多进程 拆分受体结构复合物为单独受体与配体

    Parameter
    ----------
    receptor_list : list
        受体信息列表(PDBID, 配体ID)
    '''
    cwd = get_project_dir()
    logger.debug('Current working directory: %s' % cwd)
    progress, task = _get_progress('Spliting structures', 'bold cyan', len(receptor_list))
    progress.start()
    progress.start_task(task)
    def _update(*arg):
        progress.update(task, advance=1)

    minimize_dir = cwd + '/minimize/'
    complex_dir = cwd + '/complex/'
    protein_dir = cwd + '/protein/'
    ligand_dir = cwd + '/ligands/'
    logger.debug('Splitting files in %s' % minimize_dir)
    logger.debug('Protein files will be saved in %s' % protein_dir)
    logger.debug('Ligand files will be saved in %s' % ligand_dir)

    def _move(files_tuple):
        '''
        移动文件到默认位置
        '''
        lig_file_mae = files_tuple[0]
        recep_file_mae = files_tuple[1]
        lig_file_pdb = lig_file_mae.split('.')[0] + '.pdb'
        recep_file_pdb = recep_file_mae.split('.')[0] + '.pdb'

        os.system('mv %s %s' % (lig_file_mae, ligand_dir))
        os.system('mv %s %s' % (lig_file_pdb, ligand_dir))
        os.system('mv %s %s' % (recep_file_mae, protein_dir))
        os.system('mv %s %s' % (recep_file_pdb, protein_dir))
        _update()
    # 暂时进入复合物存放的文件夹准备拆分复合物
    os.chdir(complex_dir)
    pool = multiprocessing.Pool(os.cpu_count(), maxtasksperchild=1)

    for pdbid, lig in receptor_list:
        complex_file_path = minimize_dir + '%s_minimized.mae' % pdbid
        pool.apply_async(split_com, (complex_file_path, lig), callback=_move, error_callback=error_handler)
    
    pool.close()
    pool.join()
    progress.stop()
    # 返回原始工作目录
    os.chdir(cwd)

def dock_in_pdbdir(pdbid, lig_file_path, grid_file_path, precision, calc_rmsd):
    '''
    进入PDB目录中作为对接位置的对接
    '''
    cwd = get_project_dir()
    dockfiles_dir = cwd + '/dockfiles/'
    pdb_dir = dockfiles_dir + pdbid + '/'

    mkdirs([pdb_dir])
    # 暂时进入dockfiles下的PDB文件夹 以计算对接并储存结果文件
    os.chdir(pdb_dir)
    dock(lig_file_path, grid_file_path, precision, calc_rmsd)
    os.chdir(cwd)

def self_dock(receptor_list:list, precision:str='SP', calc_rmsd:bool=True):
    '''
    共结晶配体自对接

    Parameter
    ----------
    receptor_list : list
        受体信息列表(PDBID, 配体ID)
    '''
    cwd = get_project_dir()
    logger.debug('Current working directory: %s' % cwd)
    progress, task = _get_progress('Self-docking', 'bold cyan', len(receptor_list))
    progress.start()
    progress.start_task(task)
    def _update(*arg):
        progress.update(task, advance=1)
    def _error_handler(error):
        _update()
        logger.error(error)

    ligand_dir = cwd + '/ligands/'
    grid_dir = cwd + '/grid/'
    logger.debug('Grid files will be used in %s' % grid_dir)
    logger.debug('Ligand files will be used in %s' % ligand_dir)

    pool = multiprocessing.Pool(os.cpu_count(), maxtasksperchild=1)

    for pdbid, lig in receptor_list:
        lig_file_path = ligand_dir + '%s-lig-%s.mae' % (pdbid, lig)
        grid_file_path = grid_dir + '%s_glide_grid_%s.zip' % (pdbid, lig)
        pool.apply_async(dock_in_pdbdir, (pdbid, lig_file_path, grid_file_path, precision, calc_rmsd), error_callback=_error_handler, callback=_update)
    
    pool.close()
    pool.join()
    progress.stop()

def multi_dock(mapping, precision:str='SP'):
    '''
    集合式对接核心
    
    Parameter
    ----------
    mapping : tuple|list
        映射关系元组|列表
        元组|列表中的每一个元素组成应该为(PDBID, 共结晶配体ID, 外部配体名)
    precision : str
        设定对接精度 默认SP
    '''
    cwd = get_project_dir()
    ligand_dir = cwd + '/ligands/'
    grid_dir = cwd + '/grid/'
    logger.debug('Current working directory: %s' % cwd)
    logger.debug('Grid files will be used in %s' % grid_dir)
    logger.debug('Ligand files will be used in %s' % ligand_dir)
    logger.debug('Number of all jobs: %s' % len(mapping))

    progress, task = _get_progress('Ensemble docking', 'bold cyan', len(mapping))
    progress.start()
    progress.start_task(task)
    def _update(*arg):
        progress.update(task, advance=1)
        arg.result()
    '''
    def _error_handler(error):
        _update()
        logger.error(error)
    '''

    cpu = os.cpu_count()
    logger.debug('Using Number of CPU: %s' % cpu)
    
    # ProcessPoolExecutor 需要python >= 3.3
    with ProcessPoolExecutor(cpu) as pool:
        for pdbid, self_lig, ex_lig in mapping:
            grid_file_path = grid_dir + '%s_glide_grid_%s.zip' % (pdbid, self_lig)
            lig_file_path = ligand_dir + '%s.mae' % ex_lig
            '''
            if self_lig == ex_lig:
                lig_file_path = ligand_dir + '%s_lig_%s.mae' % (pdbid, self_lig)
            '''
            future = pool.submit(dock_in_pdbdir, pdbid, lig_file_path, grid_file_path, precision, False)
            future.add_done_callback(_update)

    # Esemble dock中 multiprocessing.Pool存在重大BUG 已弃用
    # 当ligands数量较多时 父进程将卡死于pool.join()
    '''
    pool = multiprocessing.Pool(cpu, maxtasksperchild=1)

    for pdbid, self_lig, ex_lig in mapping:
        # self_lig: 该结晶自身的共结晶配体
        # ex_lig: 单次对接要用的配体
        # 当self_lig == ex_lig时即共结晶配体自对接
        grid_file_path = grid_dir + '%s_glide_grid_%s.zip' % (pdbid, self_lig)
        lig_file_path = ligand_dir + '%s.mae' % ex_lig
        if self_lig == ex_lig:
            lig_file_path = ligand_dir + '%s_lig_%s.mae' % (pdbid, self_lig)
        pool.apply_async(dock_in_pdbdir, (pdbid, lig_file_path, grid_file_path, precision, False), error_callback=_error_handler, callback=_update)

    pool.close()
    pool.join()
    # progress.update(completed = len(mapping))
    '''
    progress.stop()


def cal_mmgbsa_in_pdbdir(pdbid:str, dock_file_path:str):
    '''
    进入PDB目录并以之作为结果目录的MMGBSA能量计算
    '''
    cwd = get_project_dir()
    dockfiles_dir = cwd + '/dockfiles/'
    pdb_dir = dockfiles_dir + pdbid + '/'
    mkdirs([pdb_dir])
    # 暂时进入dockfiles下的PDB文件夹 以计算对接并储存结果文件
    os.chdir(pdb_dir)
    cal_mmgbsa(dock_file_path)
    os.chdir(cwd)

def multi_cal_mmgbsa(mapping, precision:str='SP'):
    '''
    多进程 计算多个结构的MMGBSA结合能

    Parameter
    ----------
    mapping : tuple|list
        映射关系元组|列表
        元组|列表中的每一个元素组成应该为(PDBID, 共结晶配体ID, 外部配体名)
    precision : str
        对接精度 用以定位文件
    '''    
    progress, task = _get_progress('MMGBSA Calculating', 'bold cyan', len(mapping))
    progress.start()
    progress.start_task(task)
    def _update(*arg):
        progress.update(task, advance=1)

    pool = multiprocessing.Pool(os.cpu_count(), maxtasksperchild=1)
    for pdbid, self_lig, ex_lig in mapping:
        dock_file_path = '%s_%s_glide_dock_%s_%s.maegz' % (pdbid, self_lig, ex_lig, precision)
        pool.apply_async(cal_mmgbsa_in_pdbdir, (pdbid, dock_file_path), error_callback=error_handler, callback=_update)

    pool.close()
    pool.join()
    progress.stop()

def error_handler(error):
    '''
    Error 信息显示
    '''
    logger.error(error)
    raise RuntimeError

