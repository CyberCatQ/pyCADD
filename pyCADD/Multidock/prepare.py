import os
import threading
import logging
logger = logging.getLogger('pyCADD.Multidock.prepare')

from schrodinger import structure as struc
from pyCADD.utils.check import check_file
from pyCADD.utils.getinfo import get_project_dir
from pyCADD.Dock.prepare import getpdb
from pyCADD.utils.tool import _get_progress

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
    logger.debug('Prepare to split structure file %s to %s' % (maefile, dirname))

    
    ligand_list = []
    label_list = []
    ligand_strucs = struc.StructureReader(maefile)
    st_list = [st for st in ligand_strucs]                          # 结构读入内存

    progress, task = _get_progress('Reading ligands', 'bold cyan', len(st_list))
    progress.start()
    progress.start_task(task)

    for st in st_list:
        st_name = st.property['s_m_title']                          # 分子名称
        st_activity = ''                                            # 活性标签

        # 可能的活性标签名称 bool | string
        label_name = ('b_user_Activity', 's_user_Activity', 'b_user_activity', 's_user_activity')
        for label in label_name:
            try:
                st_activity = st.property[label]
            except KeyError:
                continue
        
        # 判断是否是同名的立体异构化合物
        i = 2
        while True:
            if st_name in ligand_list:
                logger.debug('One homonymic ligand detected: %s' % st_name)
                st_name = st.property['s_m_title'] + '-%s' % str(i)
                i += 1
            else:
                break
        ligand_list.append(st_name)
        label_list.append((st_name, st_activity))
        st_maefile = dirname + st_name + '.mae'
        if not check_file(st_maefile):
            st.write(st_maefile)                                    # 写入单个mae文件
        
        progress.update(task, advance=1)

    # 创建标签文件
    with open(dirname + 'label.csv', 'w') as f:
        f.write('Ligand,activity\n')
        for st_name, st_activity in label_list:
            f.write('%s,%s\n' % (st_name, st_activity))

    progress.stop()
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
    logger.debug('PDB list: %s' % pdblist)
    logger.debug('Downloading PDB files to %s\n' % dirname)

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


