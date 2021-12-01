import logging
import multiprocessing
import os

from pyCADD.Dock.core import (dock, grid_generate, keep_chain, minimize,
                              split_com)
from pyCADD.Multidock.prepare import minimize_prepare, split_ligand
from pyCADD.utils.getinfo import get_pdbfile_path_list, get_project_dir
from pyCADD.utils.tool import mkdirs

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

    with open(list_file_path, 'r') as f:
        # 读取PDB列表中的每一行PDBID与配体名称信息
        pdbs_withlig = f.read().splitlines()        

    for i in pdbs_withlig:                                      # 按逗号分割解析列表中的PDB ID与配体名称

        pdb = i.split(',')[0].strip().upper()                   # PDB ID
        lig = i.split(',')[1].strip().upper()                   # 配体名称
        receptor_list.append((pdb, lig))                        # 将每一对作为元组储存至列表receptor_list
        
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
        for ligand in ligand_list:
            tup += ((receptor, lig, ligand),)
    
    return tup
    

def multi_minimize(pdblist:list):
    '''
    多进程 处理多个受体结构的优化
    
    Parameter
    ---------
    pdblist : list
        受体PDB ID列表
    '''
    minimize_prepare(pdblist)
    
    cwd = get_project_dir()
    # minimized文件存放目录
    minimize_dir = cwd + '/minimize/'
    pdbfiles = get_pdbfile_path_list(pdblist)

    # 暂时进入minimize文件存放目录 以优化晶体并保存结构文件于此处
    os.chdir(minimize_dir)
    # 最大进程数为CPU核心数量 1:1
    pool = multiprocessing.Pool(os.cpu_count(), maxtasksperchild=1)

    for pdbfile in pdbfiles:
        pool.apply_async(minimize, (pdbfile,), error_callback=error_handler, callback=success_handler)

    pool.close()
    pool.join()

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
    grid_dir = cwd + '/grid/'
    minimize_dir = cwd + '/minimize/'
    
    # 暂时进入grid文件存放目录 计算格点文件并保存结构于此处
    os.chdir(grid_dir)
    pool = multiprocessing.Pool(os.cpu_count(), maxtasksperchild=1)

    for pdbid, lig in receptor_list:
        st_file_path = minimize_dir + '%s_minimized.mae' % pdbid
        # 默认格点大小20A
        pool.apply_async(grid_generate, (pdbid, lig, st_file_path, 20), error_callback=error_handler, callback=success_handler)
    
    pool.close()
    pool.join()

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
    minimize_dir = cwd + '/minimize/'
    complex_dir = cwd + '/complex/'
    protein_dir = cwd + '/protein/'
    ligand_dir = cwd + '/ligands/'

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
        
    # 暂时进入复合物存放的文件夹准备拆分复合物
    os.chdir(complex_dir)
    pool = multiprocessing.Pool(os.cpu_count(), maxtasksperchild=1)

    for pdbid, lig in receptor_list:
        complex_file_path = minimize_dir + '%s_minimized.mae' % pdbid
        pool.apply_async(split_com, (complex_file_path, lig), callback=_move, error_callback=error_handler)
    
    pool.close()
    pool.join()

    # 返回原始工作目录
    os.chdir(cwd)

def dock_in_pdbdir(pdbid, lig_file_path, grid_file_path, precision, calc_rmsd):
    '''
    指定对接位置在PDB目录中的对接
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
    ligand_dir = cwd + '/ligands/'
    grid_dir = cwd + '/grid/'
    pool = multiprocessing.Pool(os.cpu_count(), maxtasksperchild=1)

    for pdbid, lig in receptor_list:
        lig_file_path = ligand_dir + '%s_lig_%s.mae' % (pdbid, lig)
        grid_file_path = grid_dir + '%s_glide_grid_%s.zip' % (pdbid, lig)
        pool.apply_async(dock_in_pdbdir, (pdbid, lig_file_path, grid_file_path, precision, calc_rmsd), error_callback=error_handler)
    
    pool.close()
    pool.join()


def multi_dock():
    '''
    集合式对接
    '''
    pass

def multi_cal_mmgbsa():
    '''
    多进程 计算多个结构的MMGBSA结合能
    '''
    pass

def read_result():
    '''
    读取对接结果文件 提取信息
    '''

def classify():
    '''
    根据信息标签 分类对接结果数据
    '''
    pass

def save_data():
    '''
    保存已分类数据
    '''
    pass


def error_handler(error):
    '''
    Error 信息显示
    '''
    logger.exception(error)

def success_handler(info):
    '''
    任务执行完成信息
    '''
    logger.debug('%s saved.' % info)
