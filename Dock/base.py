from pyCADD.utils.getinfo import downloadPDB, get_ligname, get_pdbid, catch_lig
from pyCADD.Dock import core
import os

class Dock:

    '''
    Script Core For Schrodinger Suite Analysis 
    Version 1.12

    Author: YH. W
    Last Update: 2021/11/26
        
    Parameters
    ----------
    pdbid : str
        PDB ID
    ligname : str
        配体文件PATH

    '''
    def __init__(self) -> None:

        self.pdbid = get_pdbid()
        downloadPDB(self.pdbid)
        self.pdbfile = self.pdbid + '.pdb'           # PDB结构文件
        self.ligname = get_ligname()

        self.minimized_file = ''                     # Minimized化完成的文件名
        self.grid_file = ''                          # 格点文件名
        self.lig_file = ''          # 内源配体文件名
        self.recep_file = ''        # 受体文件名
        self.dock_file = ''         # 对接结果文件名
        self.mmgbsa_file = ''       # 结合能计算结果文件名
        self.sitemap_file = ''      # 结合口袋体积计算结果文件名
    
    def __if_keep_single_chain(self) -> str:
        '''
        询问是否保留单链

        Return
        ----------
        str
            处理完成的PDB结构文件名
        '''
        lis = catch_lig(self.pdbfile)
        if len(lis) > 1:
            print('\n')
            os.system('cat %s | grep -w -E ^HET' % self.pdbfile)
            print('There are multiple ligand small molecules. \nDo you need to keep a single chain？(Y/N)')
            _flag = input().strip().upper() # 是否保留单链的标志

            if _flag == 'Y':
                chain = input('Enter the Chain Code:').strip().upper()
                pdbfile = core.keep_chain(self.pdbfile, chain)
                self.pdbfile = pdbfile  # [属性修改]修改pdbfile为单链文件
                return pdbfile
            else:
                print('\nChoosed Original Crystal.\n')
                return self.pdbfile