import re
import os

def check_file(file_path:str) -> bool:
    '''
    检查文件是否存在
    存在返回True 否则返回False

    Parameter
    ----------
    file_path : str
        文件PATH
    
    Return
    ----------
    '''
    return os.path.exists(file_path)

def checkpdb(pdb: str):
    '''
    检查PDB ID合法性

    Return
    -------
    bool
        合法返回True 否则返回False
    '''

    match = re.fullmatch(r'^\d[0-9a-zA-Z]{3,}$', pdb)
    if match:
        return True
    else:
        return False


def check_ligname(ligname: str):
    '''
    检查配体名称合法性
    Parameter
    ----------
    ligname : str
        配体名称

    Return
    ----------
    match[str] | False
        合法返回Match对象 否则返回False
    '''

    match = re.search('[0-9A-Z]{2,3}$', ligname)
    if match:
        return match
    else:
        return False

def check_chain(pdbfile:str) -> str:
    '''
    检查链数并询问是否保留单链

    Parameter
    ---------
    pdbfile : str
        PDB文件PATH

    Return
    ----------
    str
        处理完成的PDB结构文件名
    '''
    from pyCADD.utils.getinfo import catch_lig
    from pyCADD.Dock.core import keep_chain

    lis = catch_lig(pdbfile)
    if len(lis) > 1:
        print('\n')
        os.system('cat %s | grep -w -E ^HET' % pdbfile)
        print('There are multiple ligand small molecules. \nDo you need to keep a single chain？(Y/N)')
        _flag = input().strip().upper() # 是否保留单链的标志

        if _flag == 'Y':
            chain = input('Enter the Chain Code:').strip().upper()
            pdbfile = keep_chain(pdbfile, chain)

            return pdbfile
        else:
            print('\nChoosed Original Crystal.\n')
            return pdbfile
    
    else:
        return pdbfile

def check_lig_num(maefile:str) -> int:
    '''
    检查maestro文件中含有多少个小分子化合物

    Parameter
    ---------
    maefile : str
        maestro文件PATH

    Return
    ---------
    int
        化合物数量
    '''
    from schrodinger import structure as struc
    ligand_strucs = struc.StructureReader(maefile)
    
    count = 0
    for st in ligand_strucs:
        count += 1
    
    return count


