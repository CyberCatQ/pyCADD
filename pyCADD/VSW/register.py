import logging
import os

logger = logging.getLogger('pyCADD.VSW.register')
vsw_dir = os.path.dirname(os.path.abspath(__file__)) + '/'

from pyCADD.utils.tool import Myconfig

gene_config_file = vsw_dir + 'genelist.ini'
database_config_file = vsw_dir + 'database.ini'

config_gene = Myconfig()
config_database = Myconfig()

config_gene.read(gene_config_file)
config_database.read(database_config_file)

def reg_gene(gene:str, family:str='GENE', path:str=''):
    '''
    注册基因并添加到配置文件中

    Parameters
    ----------
    gene : str
        注册基因名
    family : str
        基因家族标签名(可选)
    path : str
        基因PDBID 列表文件路径
    '''

    if config_gene.has_section(family):
        logger.debug('Section %s current gene list: %s' % (family, config_gene.options(family)))
    else: 
        logger.debug('Create new section: %s' % family)

    # 将其拷贝至程序目录进行注册
    os.system('cp %s %s' % (path, vsw_dir))
    value = vsw_dir + os.path.basename(path)

    config_gene.set(family, gene, value)
    config_gene.write(open(gene_config_file, 'w'))
    logger.debug('Section %s current gene list: %s' % (family, config_gene.options(family)))
    logger.info('Gene %s has been registed to %s.' % (gene, family))

def reg_database(database:str, label:str='DATABASE', path:str=''):
    '''
    注册基因并添加到配置文件中

    Parameters
    ----------
    database : str
        注册数据库名称
    family : str
        数据库标签
    path : str
        注册数据库路径
    '''

    config_database.set(label, database, path)
    config_database.write(open(database_config_file, 'w'))
    logger.info('Database %s has been registed to %s.' % (database, label))
    
def del_gene(gene:str, family:str='GENE'):
    '''
    从注册配置文件中删除指定基因

    Parameters
    ----------
    gene : str
        注册基因名
    family : str
        基因家族标签名(可选)
    '''
    if not config_gene.has_option(family, gene):
        logger.error('%s is not found in %s' % (gene, family))
    else:
        config_gene.remove_option(family, gene)
        logger.info('%s in %s has been removed.' % (gene, family))
    config_gene.write(open(gene_config_file, 'w'))
    
def del_database(database:str, label:str='DATABASE'):
    '''
    从注册配置文件中删除指定数据库

    Parameters
    ----------
    database : str
        注册数据库名
    label : str
        数据库分类标签
    '''
    if not config_database.has_option(label, database):
        logger.error('%s is not found in %s' % (database, label))
    else:
        config_database.remove_option(label, database)
        logger.info('%s in %s has been removed.' % (database, label))
    config_database.write(open(database_config_file, 'w'))
