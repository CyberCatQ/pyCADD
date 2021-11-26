import re

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
