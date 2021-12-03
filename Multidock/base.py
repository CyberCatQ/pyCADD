import os
import logging

from pyCADD.utils.tool import generate_logfile_name, mkdirs
from pyCADD.Multidock import core
from pyCADD.utils.getinfo import get_pdblist_from_recplist, get_project_dir

base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# log_dir = base_dir + '/logs'

#  配置log
logger = logging.getLogger('pyCADD')
logger.setLevel(level = logging.INFO)
file_fmt = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
console_fmt = logging.Formatter('%(levelname)s - %(message)s')

global logfile
logfile = generate_logfile_name()
filehandler = logging.FileHandler(logfile, 'a')
filehandler.setLevel(logging.INFO)
filehandler.setFormatter(file_fmt)
consolehandler = logging.StreamHandler()
consolehandler.setLevel(logging.INFO)
consolehandler.setFormatter(console_fmt)

logger.addHandler(filehandler)
logger.addHandler(consolehandler)


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

    @property
    def logfile(self):
        '''
        获取当前任务log文件PATH
        '''
        global logfile
        return  logfile

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
   
    def map(self):
        '''
        建立映射关系
        '''
        self.mapping = core.map(self.receptor_list, self.ligand_list)
        logger.info('Create mapping complete.')

    def get_pdblist(self):
        '''
        获取PDB ID列表
        '''
        return self.pdblist
    
    def optimize(self):
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

    def ensemble_dock(self, precision:str='SP'):
        '''
        启动集合式对接
        '''
        logger.info('Prepare to run ensemble docking')
        core.multi_dock(self.mapping, precision)
        logger.info('Ensemble Docking Done.')
    
    def calc_mmgbsa(self, precision:str='SP'):
        '''
        计算MMGBSA结合能
        '''
        logger.info('Prepare to calculate MM-GB/SA Binding Energy')
        logger.warning('Attention: It may cost extremely much time')
        core.multi_cal_mmgbsa(self.mapping, precision)
        logger.info('MM-GB/SA binding energy calculate Done')
    

