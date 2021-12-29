import os
from schrodinger import structure as struc
from pyCADD.utils.getinfo import downloadPDB as get_pdb
import logging
logger = logging.getLogger('pyCADD.Dock.prepare')

def load_st(st_file: str) -> object:
    '''
    读取结构(如有多个结构仅返回第一个)

    Parameters
    ----------
    st_file : str
        需要读取结构的文件path 需要Schrodinger支持的文件格式

    Return
    ---------
    object
        结构对象
    '''
    logger.debug('Loading structure file: %s' % st_file)
    return next(struc.StructureReader(st_file))

def convert_format(file_path:str, suffix:str) -> str:
    '''
    mae/pdb格式转换为pdb/mae格式

    Parameters
    ----------
    file_path : str
        需要转换的文件PATH
    suffix : str
        格式后缀名(pdb|mae)

    Return
    ----------
    str
        转换后的文件名
    '''
    st = load_st(file_path)
    file = os.path.basename(file_path)
        
    if suffix == 'pdb':
        convert_file = file.split('.')[0] + '.pdb'
    elif suffix == 'mae':
        convert_file = file.split('.')[0] + '.mae'

    st.write(convert_file)
    return convert_file

def getpdb(pdbid:str) -> str:
    '''
    下载PDB文件
    '''
    logger.info('Downloading PDB file: %s' % pdbid + '.pdb')
    return get_pdb(pdbid)
