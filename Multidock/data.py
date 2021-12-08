import multiprocessing
import os
import logging

import pandas as pd
from pyCADD.Dock.data import extra_data
from pyCADD.utils.getinfo import get_project_dir
from pyCADD.utils.tool import _get_progress

logger = logging.getLogger('pyCADD.Multidock.data')

def read_reuslt(pdbid:str, result_file_path:str):
    '''
    进入PDB目录并读取结果文件数据
    '''
    cwd = get_project_dir()
    dockfiles_dir = cwd + '/dockfiles/'
    pdb_dir = dockfiles_dir + pdbid + '/'

    # 暂时进入PDB文件夹查找结果文件
    os.chdir(pdb_dir)
    data_dic = extra_data(result_file_path)
    os.chdir(cwd)

    return data_dic

def multi_read_result(mapping, precision:str='SP', mmgbsa:bool=False):
    '''
    多进程 读取对接结果文件 提取信息

    Parameter
    ---------
    mapping : tuple|list
        映射关系元组|列表
        元组|列表中的每一个元素组成应该为(PDBID, 共结晶配体ID, 外部配体名)
    precision : str
        对接精度
    mmgbsa : bool
        是否提取MMGBSA结合能计算结果
    '''
    progress, task = _get_progress('Extracting data', 'bold cyan', len(mapping))
    progress.start()
    progress.start_task(task)
    data_list = []

    def _update(*arg):
        progress.update(task, advance=1)
    def _append_to_lis(dic):
        _update()
        if dic:
            data_list.append(dic)

    pool = multiprocessing.Pool(os.cpu_count(), maxtasksperchild=1)
    with_mmgbsa = ''
    if mmgbsa:
        with_mmgbsa = '_mmgbsa'
        
    for pdbid, self_lig, ex_lig in mapping:
        dock_file_path = '%s_%s_glide_dock_%s_%s%s.maegz' % (pdbid, self_lig, ex_lig, precision, with_mmgbsa)
        pool.apply_async(read_reuslt, (pdbid, dock_file_path,), callback=_append_to_lis)

    pool.close()
    pool.join()
    progress.stop()

    return data_list

def save_data(data_list:list, pdb_list:list):
    '''
    按照PDB 分类保存数据
    '''
    cwd = get_project_dir()
    result_dir = cwd + '/results/'
    logger.debug('Results files will be saved in: %s' % result_dir)
    os.chdir(result_dir)

    data = pd.DataFrame(data_list)
    data.to_csv('TOTAL.csv')
    group = data.groupby('PDB')
    for pdb in pdb_list:
        group.get_group(pdb).to_csv('%s_DOCK_FINAL_RESULTS.csv' % pdb, index=False)
        logger.debug('PDB %s data saved.' % pdb)
    

