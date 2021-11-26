import os
from pyCADD.utils.check import checkpdb, check_ligname
root_path = os.path.abspath(os.path.dirname(__file__)).split('pyCADD')[0]  # 项目总路径
base_path = root_path + 'pyCADD/'                          # pyCADD程序路径
base_url = 'https://files.rcsb.org/download/'

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

    print('Downloading %s ...' % pdbid)
    url = base_url + pdbid + '.pdb'
    os.system('wget -q %s' % url)

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
            pdbid = str(input('\nTo get PDB ID automatically, please change the name of crystal folder to PDBID\n Input PDB ID:')).strip().upper()
            if checkpdb(pdbid):
                return pdbid
            else:
                print('请输入正确的PDB ID!')

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

def get_ligname(pdbid=None) -> str:
    '''
    尝试自动获取配体名 如有多个配体则获取用户输入的配体名并检查合法性

    Return
    ----------
    str
        自动识别或手动输入的配体名称

    '''

    if not pdbid:
        pdbid = get_pdbid()
    
    pdbfile = pdbid + '.pdb'

    if not os.path.exists(pdbfile):  # 下载PDB文件

        print('%s.pdb is not found. Would you like to download it? (Y/N)' % pdbid)
        if input().upper().strip() == 'Y':
            downloadPDB(pdbid)
        else:
            print('PDB File is not Found. Exit.')

    # 抓取pdb原始结构文件(不可是已处理过的结构)中单一entry关于小分子的描述
    lis = catch_lig(pdbfile)  
    lig = []
    if len(lis) == 0:
        raise RuntimeError('%s is an Apo Crystal.' % pdbid)

    for i in lis:
        passed = check_ligname(i.strip())  # 匹配得到配体候选列表
        if passed:
            lig_name = passed.group()
            if lig_name != 'HOH' and not lig_name in lig:  # 排除配体候选中的水分子与重复分子
                lig.append(lig_name)

    if len(lig) == 1:
        ligname = str(lig[0])
        return ligname
    else:
        print('Crystal %s has more than one ligand:' % pdbid, ''.join(str(x)+' ' for x in lig), end='\n')
        while True:
            ligname = input('Please specify ligand name:').strip().upper()
            if check_ligname(ligname) and ligname in lig:
                return ligname
            else:
                print('Wrong ligand name, please try again.')
    