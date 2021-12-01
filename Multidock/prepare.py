import os
import threading

from schrodinger import structure as struc
from pyCADD.utils.check import check_file
from pyCADD.utils.getinfo import get_project_dir
from pyCADD.Dock.prepare import getpdb

def split_ligand(maefile:str, dirname:str='./') -> list:
    '''
    将单个mae文件中包含的所有小分子 每一个拆分为独立mae文件

    Parameter
    ----------
    maefile : str
        要拆解的mae文件PATH
    dirname : str
        拆解的小分子文件存放目录PATH
    
    Return
    ----------
    list
        小分子名称列表
    '''
    dirname = os.path.abspath(dirname) + '/'
    ligand_list = []
    ligand_strucs = struc.StructureReader(maefile)

    for st in ligand_strucs:
        st_name = st.property['s_m_title']                          # 分子名称
        st_maefile = dirname + st_name + '.mae'
        if not check_file(st_maefile):
            st.write(st_maefile)                        # 写入单个mae文件
        ligand_list.append(st_name)
    
    return ligand_list

def download_pdblist(pdblist:list, dirname:str):
    '''
    多线程下载PDB列表中的文件
    Parameters
    ----------
    pdblist : list
        要下载的PDB列表
    dirname : str
        保存PDB文件的目录PATH
    '''

    dirname = os.path.abspath(dirname) + '/'
    pdblist_to_download = []
    cwd = os.getcwd()

    for pdb in pdblist:
        pdbfile = dirname + pdb + '.pdb'
        if not check_file(pdbfile):
            pdblist_to_download.append(pdb)

    # 暂时进入PDB存放目录 以下载文件到此处
    os.chdir(dirname)
    threads = []
    for pdb in pdblist_to_download:
        t = threading.Thread(target=getpdb, args=(pdb,))
        threads.append(t)
    
    for t in threads:
        t.start()
    for t in threads:
        t.join() 

    # 返回默认工作目录
    os.chdir(cwd)

def minimize_prepare(pdblist:list):
    '''
    多线程结构优化准备工作

    Parameter
    ----------
    pdb_list : list
        受体PDBID列表
    
    '''
    cwd = get_project_dir()
    pdb_dir = cwd + '/pdb/'
    download_pdblist(pdblist, pdb_dir)


