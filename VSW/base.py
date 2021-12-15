import logging
import os

from pyCADD.Multidock.base import Multidock
from pyCADD.utils.check import check_file
from pyCADD.utils.getinfo import get_project_dir
from pyCADD.utils.tool import mkdirs
from pyCADD.VSW import core
from rich.prompt import Prompt

logger = logging.getLogger('pyCADD.VSW.base')

vsw_dir = os.path.dirname(os.path.abspath(__file__)) + '/'
gene_config_path = vsw_dir + 'genelist.ini'
database_config_path = vsw_dir + 'database.ini'

class VSW(Multidock):
        
    '''
    Python Script for Virtual Screening Workflow.
    '''

    def __init__(self) -> None:
        self.pdblist = []               # 由基因索取到的PDB ID列表 由元组(PDBID, Ligand)组成
        self.grid_num = 1               # 总计需要进行VSW的PDB数量(创建的pipeline数量)
        self.gene_config = {}           # 受体配置信息
        self.database_config = {}       # 化合物库配置信息
        mkdirs(self.required_dir)
        self.read_gene()
        self.read_databse()

    @property
    def vsw_dir(self) -> str:
        return self.project_dir + '/vsw/'

    @property
    def required_dir(self):
        return [
            self.complex_dir, 
            self.dockfiles_dir, 
            self.grid_dir, 
            self.ligands_dir, 
            self.minimize_dir, 
            self.pdb_dir, 
            self.protein_dir, 
            self.log_dir, 
            self.result_dir,
            self.vsw_dir
            ]

    @property
    def genelist(self) -> list:
        '''
        基因成员列表
        '''
        return list(self.gene_config.keys())
    
    @property
    def database_list(self) -> list:
        '''
        化合物库成员列表
        '''
        if not self.database_config:
            raise ValueError('No database config loaded.')
        return list(self.database_config.keys())

    def read_gene(self) -> str:
        '''
        读取受体信息
        '''
        logger.info('Reading gene info from %s' % gene_config_path)
        self.gene_config = core.read_gene_config(gene_config_path)
        logger.info('Loaded gene info.')

    def read_databse(self) -> str:
        '''
        读取化合物库信息
        '''
        logger.info('Reading database info from %s' % database_config_path)
        self.database_config = core.read_database_config(database_config_path)
        logger.info('Database info loaded.')

    def get_gene(self) -> str:
        '''
        打印当前基因信息并获取用户指定的VSW基因
        '''
        logger.info('Current Gene List: %s' % self.genelist)
        self.gene = Prompt.ask('Enter the required gene', choices=self.genelist, show_choices=False)
        logger.info('Current Gene: %s' % self.gene)
        return self.gene

    def get_receptor_list(self) -> None:
        '''
        读取基因对应的PDB列表文件
        '''
        if not self.gene:
            self.get_gene()

        pdblist_file = self.vsw_dir + self.gene + '.txt'
        if not check_file(pdblist_file):
            pdblist_file = input('Enter the path of %s PDB list: ' % self.gene).strip()
            if not check_file(pdblist_file):
                raise(FileNotFoundError('%s not found.' % pdblist_file))

        self.read_receptor(pdblist_file)
        logger.info('Current Gene PDB List: %s' % self.pdblist)
        logger.info('PDB list length: %s' % str(len(self.pdblist)))
    
    def select_database(self) -> None:
        '''
        打印当前化合物库信息并获取用户指定的化合物库名称
        '''
        logger.info('Current compounds library list: %s' % self.database_list)
        self.database = Prompt.ask('Enter the required library', choices=self.database_list, show_choices=False)
        logger.info('Current database: %s' % self.database)
        logger.info('Current database file path: %s' % self._get_database_path(self.database))
        return self.database
    
    def _get_database_path(self, name) -> None:
        '''
        返回化合物库路径
        '''
        return self.database_config[name]

    def generate_input_file(self) -> None:
        '''
        生成VSW输入文件
        '''
        jobname = self.gene + '_' + 'VSW'
        core.gen_input_file(self.receptor_list, self._get_database_path(self.database), jobname)
        logger.info('VSW input file %s generated.' % (jobname + '.inp'))

    
