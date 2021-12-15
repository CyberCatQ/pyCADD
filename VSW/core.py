from pyCADD.Dock.core import launch
from pyCADD.utils.check import check_file
from pyCADD.utils.getinfo import get_config

import logging

logger = logging.getLogger('pyCADD.VSW.core')

def read_gene_config(gene_config_path):
    '''
    从文件中读取可供筛选的受体基因

    Parameter
    ----------
    gene_config_path : str
        受体信息配置文件路径

    Return
    ----------
    dict
        受体配置信息
    '''

    return dict(get_config(gene_config_path)['GENE'])

def read_database_config(database_config_path):
    '''
    从文件中读取化合物库路径信息

    Parameter
    ----------
    database_config_path : str
        化合物库路径信息文件路径
    
    Return
    ----------
    dict
        化合物库路径配置信息
    '''

    return dict(get_config(database_config_path)['DATABASE'])



