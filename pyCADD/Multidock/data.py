import multiprocessing
import os
import logging
import re

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

def save_data(data_list:list, pdblist:list, additional_col:list=[]):
    '''
    按照PDB 分类保存数据
    Parameters
    ----------
    data_list : list
        由数据字典组成的列表
    pdblist : list
        PDB ID列表(用以聚类结果)
    additional_col : list
        额外需要提取的数据列名称(maestro原始列名 e.g: r_i_glide_gscore)
    '''
    cwd = get_project_dir()
    result_dir = cwd + '/results/'
    logger.debug('Results files will be saved in: %s' % result_dir)
    os.chdir(result_dir)

    prop = ['PDB', 'Ligand', 'Original', 'Docking_Score', 'MMGBSA_dG_Bind', 'rmsd', 'precision', 'Site_Score', 'Volume','ligand_efficiency', 
            'XP_Hbond', 'rotatable_bonds', 'ecoul', 'evdw', 'emodel', 'energy', 'einternal', 'activity' ,'lipo', 'hbond', 'metal', 'rewards', 'erotb', 'esite']
    if additional_col:
        for col in additional_col:
            prop.append(col)
            
    data = pd.DataFrame(data_list, columns=prop)
    data.to_csv('TOTAL.csv')
    group = data.groupby('PDB')
    for pdb in pdblist:
        group.get_group(pdb).to_csv('%s_DOCK_FINAL_RESULTS.csv' % pdb, index=False)
        logger.debug('PDB %s data saved.' % pdb)

    os.chdir(cwd)
    

def merge_data(pdblist):
    '''
    抽取所有配体在全部PDB中的对接分数 并合并为矩阵

    Parameter
    ----------
    pdblist : list
        PDB ID列表
    '''
    cwd = get_project_dir()
    result_dir = cwd + '/results/'
    logger.debug('Results files will be used in: %s' % result_dir)
    os.chdir(result_dir)

    data_list = []
    for pdb in pdblist:
        resultfile = '%s_DOCK_FINAL_RESULTS.csv' % pdb
        data = pd.read_csv(resultfile, index_col='Ligand')
        ds_data = data[['Docking_Score']]
        data_list.append(ds_data.rename(columns={'Docking_Score' : pdb}))

    original_ligand_label = pd.read_csv(cwd + '/ligands/original.csv', index_col='Ligand')
    label_data = pd.read_csv(cwd + '/ligands/label.csv', index_col='Ligand')
    label_data = pd.concat([original_ligand_label, label_data])

    final = data_list[0].join(data_list[1:], how='outer')
    final = final.join(label_data, how='left')
    final.to_csv('matrix.csv', index=True)

    _IDs = {}
    _ref_ligands = {}
    for index, row in final[final['activity'] == 'origin'].iterrows():
        _IDs[index[:4]] = row[index[:4]]
        _ref_ligands[index[:4]] = re.search(r'(?<=-lig-)[a-zA-Z0-9]+', index).group()
    
    pd.DataFrame([_IDs, _ref_ligands], index=['Reference', 'Original ligand']).to_csv('reference.csv')

    os.chdir(cwd)



