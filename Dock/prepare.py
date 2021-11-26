import os
from schrodinger import structure as struc

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
        st.write(file.split('.')[0] + '.pdb')
        return file.split('.')[0] + '.pdb'
    elif suffix == 'mae':
        st.write(file.split('.')[0] + '.mae')
        return file.split('.')[0] + '.mae'

