import os

from schrodinger import structure as struc


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
        st_name = st.property['s_m_title']              # 分子名称
        st.write(dirname + st_name + '.mae')                      # 写入单个mae文件
        ligand_list.append(st_name)
    
    return ligand_list

