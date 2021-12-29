import os
import logging
import requests

from pyCADD.utils.check import checkpdb, check_ligname
from pyCADD.utils.tool import Myconfig

root_path = os.getcwd()

base_url = 'https://files.rcsb.org/download/'
logger = logging.getLogger('pyCADD.getinfo')

def downloadPDB(pdbid) -> None:
    '''
    从RCSB服务器下载PDB文件

    Parameter
    ----------
    pdbid : str
        PDB ID字符串
    '''
    pdbfile = pdbid + '.pdb'
    if os.path.exists(pdbfile):
        return

    logger.debug('Downloading %s ...' % pdbid)
    url = base_url + pdbid + '.pdb'
    response = requests.get(url)
    pdb_data = response.text
    with open(pdbid + '.pdb', 'w') as f:
        f.write(pdb_data)
    
    logger.debug('%s.pdb downloaded.' % pdbid)

def get_pdbid() -> str:
    '''
    获取用户输入的PDBID并检查合法性

    Return
    ----------
    str
        PDB ID字符串

    '''

    pdbid = root_path.split('/')[-2]

    if checkpdb(pdbid):        
        return pdbid
    else:
        while True:
            logger.info('To get PDB ID automatically, please change the name of crystal folder to PDBID.')
            pdbid = input('Input PDB ID:').strip().upper()
            logger.info('PDB ID: %s' % pdbid)
            if checkpdb(pdbid):
                return pdbid
            else:
                logger.warning('请输入正确的PDB ID!')

def catch_lig(pdbfile) -> list:
    '''
    从PDB文件获取配体小分子信息
    按行分割并返回列表

    Parameter
    ----------
    pdbfile : str
        PDB文件PATH

    Return
    ----------
    list
        配体小分子信息列表
    '''
    return  os.popen("cat %s | grep -w -E ^HET | awk '{print $2}'" % pdbfile).readlines()

def get_ligname(pdbid) -> str:
    '''
    尝试自动获取配体RCSB ID 如有多个配体则获取用户输入的配体名并检查合法性

    Return
    ----------
    str
        自动识别或手动输入的配体名称

    '''

    pdbfile = pdbid + '.pdb'

    # 抓取pdb原始结构文件(不可是已处理过的结构)中单一entry关于小分子的描述
    lis = catch_lig(pdbfile)  
    lig = []
    if len(lis) == 0:
        raise RuntimeError('%s is an Apo Crystal.' % pdbid)

    # 匹配得到配体候选列表
    for i in lis:
        passed = check_ligname(i.strip())  
        if passed:
            lig_name = passed.group()
            # 排除配体候选中的水分子与重复分子
            if lig_name != 'HOH' and not lig_name in lig:  
                lig.append(lig_name)

    if len(lig) == 1:
        ligname = str(lig[0])
        return ligname
    else:
        logger.info('Crystal %s has more than one ligand: %s' % (pdbid, ''.join(str(x)+' ' for x in lig)))
        while True:
            ligname = input('Please specify ligand name:').strip().upper()
            if check_ligname(ligname) and ligname in lig:
                return ligname
            else:
                logger.warning('Wrong ligand name, please try again.')
    
    
def get_ligmol_info(file:str, ligname:str) -> str:
    '''
    以ligname为KEY 查找Maestro文件中的Molecule Number

    Parameters
    ----------
    file : str
        maestro文件PATH
    ligname : str
        配体ID(RCSB ID)

    Return
    ----------
    str
        配体所在Molecule Number

    '''

    from pyCADD.Dock.prepare import load_st
    # 载入结构对象
    st = load_st(file)  

    # 获取结构中的配体所在Molecule object
    def _get_mol_obj(st):   
        residues = st.residue
        for res in residues:
            if res.pdbres.strip() == '%s' % ligname:
                molnum = res.molecule_number
                yield st.molecule[molnum]

    mol = next(_get_mol_obj(st))

    # 判断该molecule是否仅包括小分子本身(是否存在共价连接) 自动移除共价连接
    if len(mol.residue) != 1:  
        logger.info('%s in %s : A covalent bond may exist between the ligand and residue.' % (ligname, file))
        logger.info('An attempt will be made to remove the covalent bond automatically.')
        bonds = st.bond

        for bond in bonds:
            resname1 = bond.atom1.getResidue().pdbres.strip()
            resname2 = bond.atom2.getResidue().pdbres.strip()
            if resname1 == '%s' % ligname or resname2 == '%s' % ligname:
                if resname1 != resname2:
                    bond_to_del = bond
            
        if not bond_to_del:
            raise RuntimeError('Can not delete covalent bonds automatically.')
            
        st.deleteBond(bond_to_del.atom1, bond_to_del.atom2)
        st.write(file)

        mol = next(_get_mol_obj(st))

    return mol.number

def get_base_dir():
    '''
    获取pyCADD项目文件夹所在的Absolute PATH
    '''
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def get_project_dir():
    '''
    获取项目文件夹的Absolute PATH
    '''

    return os.getcwd()

def get_pdblist_from_recplist(receptor_list:list):
    '''
    获取输入的受体列表信息 剥离PDB ID信息并返回列表

    Parameter
    ----------
    receptor_list : list
        (受体PDBID, 配体ID)组成的列表
    
    Return
    ----------
    list
        PDB ID列表
    '''

    pdblist = []
    for pdb, lig in receptor_list:
        pdblist.append(pdb)
    
    return pdblist

def get_pdbfile_path_list(pdblist):
    '''
    将PDBID列表转换为对应的PDB文件PATH列表

    Parameter
    ----------
    pdblist : list
        PDB ID列表
    
    Return
    ----------
    list
        PDB文件PATH列表
    '''
    
    pdbfile_path_list = []
    pdb_dir = get_project_dir() + '/pdb/'
    for pdb in pdblist:
        pdbfile = pdb_dir + pdb + '.pdb'
        pdbfile_path_list.append(pdbfile)

    return pdbfile_path_list

def get_gridfile_path_list(receptor_list:list):
    '''
    将receptor list转换为对应的PDB grid文件PATH列表

    Parameter
    ----------
    receptor_list : list
        (受体PDBID, 配体ID)组成的列表
    
    Return
    ----------
    list
        Grid文件PATH列表
    '''

    gridfile_path_list = []
    grid_dir = get_project_dir() + '/grid/'
    for pdbid, lig in receptor_list:
        gridfile = grid_dir + '%s_glide_grid_%s.zip' % (pdbid, lig)
        gridfile_path_list.append(gridfile)
    
    return gridfile_path_list

def get_config(config_file_path):
    '''
    获取配置文件信息

    Parameter
    ---------
    config_file_path : str
        配置文件路径

    Return
    ---------
    ConfigParser
        配置信息对象
    '''

    config = Myconfig()
    config.read(config_file_path)

    return config
