import os

base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
log_dir = base_dir + '/logs'

from datetime import datetime

date = datetime.now()
year = str(date.year)
month = str(date.month)
day = str(date.day)
now = year + month.rjust(2, '0') + day.rjust(2, '0')

# 配置log
import logging

logger = logging.getLogger('pyCADD')
logger.setLevel(level = logging.INFO)
file_fmt = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
console_fmt = logging.Formatter('%(levelname)s - %(message)s')

filehandler = logging.FileHandler(log_dir + '/%s.log' % now, 'a')
filehandler.setLevel(logging.INFO)
filehandler.setFormatter(file_fmt)
consolehandler = logging.StreamHandler()
consolehandler.setLevel(logging.INFO)
consolehandler.setFormatter(console_fmt)

logger.addHandler(filehandler)
logger.addHandler(consolehandler)

from pyCADD.Multidock import core
from pyCADD.utils.getinfo import get_pdblist_from_recplist, get_project_dir
from pyCADD.utils.tool import mkdirs


class Multidock:
    '''
    Multi Mode
    '''

    def __init__(self) -> None:
        self.pdblist = []
        self.receptor_list = []
        self.ligand_list = []
        mkdirs(self.required_dir)

    # 一些必要的目录
    @property
    def project_dir(self):
        return get_project_dir()

    @property
    def dockfiles_dir(self):
        return self.project_dir + '/dockfiles'
    
    @property
    def grid_dir(self):
        return self.project_dir + '/grid'

    @property 
    def ligands_dir(self):
        return self.project_dir + '/ligands'

    @property
    def complex_dir(self):
        return self.project_dir + '/complex'

    @property
    def minimize_dir(self):
        return self.project_dir + '/minimize'
    
    @property
    def protein_dir(self):
        return self.project_dir + '/protein'
    
    @property
    def pdb_dir(self):
        return self.project_dir + '/pdb'
        
    @property
    def required_dir(self):
        return [self.complex_dir, self.dockfiles_dir, self.grid_dir, self.ligands_dir, self.minimize_dir, self.pdb_dir, self.protein_dir]

    def read_receptor(self, file_path:str):
        '''
        读入受体列表文件
        '''
        self.receptor_list = core.read_receptors(file_path)
        self.pdblist = get_pdblist_from_recplist(self.receptor_list)
        logger.info('Read receptors file: %s' % file_path)
    
    def read_ligands(self, file_path:str):
        '''
        读入小分子结构文件(mae格式)
        '''
        self.ligand_list = core.read_ligands(file_path, self.ligands_dir)
        logger.info('Loaded ligands file: %s' % file_path)
    
    def _map(self):
        '''
        建立映射关系
        '''
        self.mapping = core.map(self.receptor_list, self.ligand_list)
        logger.info('Map complete.')

    def get_pdblist(self):
        '''
        获取PDB ID列表
        '''
        return self.pdblist
    
    def minimize(self):
        '''
        优化全部晶体
        '''
        logger.info('Prepare to minimize crystals')
        core.multi_minimize(self.pdblist)
        logger.info('Minimize Done.')

    def grid_generate(self):
        '''
        生成格点文件
        '''
        logger.info('Prepare to generate grid files')
        core.multi_grid_generate(self.receptor_list)
        logger.info('Grid generate Done.')
    
    def split(self):
        '''
        拆分复合物
        '''
        logger.info('Prepare to split complex structures')
        core.multi_split(self.receptor_list)
        logger.info('Split Done.')

    def self_dock(self, precision:str='SP', calc_rmsd:bool=True):
        '''
        共结晶配体自动对接
        '''
        logger.info('Prepare to re-dock ligand')
        core.self_dock(self.receptor_list, precision, calc_rmsd)
        logger.info('Self-Dock Done.')

    
    