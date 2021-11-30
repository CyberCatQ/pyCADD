from pyCADD.Multidock.prepare import split_ligand

def read_receptors(list_file_path:str):
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

def read_ligands(maefile:str, dirname:str='./'):
    '''
    处理配体 拆分mae文件为单个化合物的结构

    Parameter
    ---------
    maefile : str
        包含所有配体小分子的单个mae文件PATH
    dirname : str
        拆分的mae储存目录PATH
    '''

    return split_ligand(maefile, dirname)


def map(receptor_list:list, ligand_list:list):
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
    

def multi_minimize():
    '''
    多进程 处理多个受体结构的优化
    '''
    pass

def multi_grid_generate():
    '''
    多进程 生成多个受体结构的格点文件
    '''
    pass

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


    