from pyCADD.utils.check import check_file
import os
import logging

from rcsbsearch import rcsb_attributes as attrs
logger = logging.getLogger('pyCADD.query.core')

def _precess_inputfile(input_file_path:str):
    '''
    解析输入文件
    Parameter
    ----------
    input_file_path : str
        输入文件路径
    
    Return
    ----------
    dict

    '''
    if not check_file(input_file_path):
        raise FileNotFoundError('File %s not found.' % input_file_path)


    gene_dic = {}
    with open(input_file_path) as f:
        gene_infos = f.read().splitlines()

    for gene_info in gene_infos:
        gene_info = gene_info.strip().split(',')
        gene = gene_info[0]
        name = gene_info[1]
        gene_dic[gene] = name

    return gene_dic

def get_input():
    '''
    读取并判断输入类型 返回基因表

    Return
    ----------
    dict
        基因 : 常用名
    '''
    input_text = input('Enter the gene name OR path of gene list: ')
    if os.path.isfile(input_text):
        logger.debug('Input a file path: %s' % input_text)
        logger.debug('Process as input file.')
        return _precess_inputfile(input_text)
    else:
        logger.debug('Input a gene-like string: %s' % input_text)
        gene, name = input_text.split(',')
        logger.debug('Gene: %s' % gene)
        logger.debug('Abbreviation: %s' % name)
        d = {}
        d[gene] = name
        return d
    
def search(gene_list=[]):
    '''
    在RCSB服务器搜索列表中的基因相关晶体PDBID 存储为文本文件

    Parameters
    ----------
    gene_list : list
        要搜索的全部基因列表
    '''


    logger.debug('\nAll Genes:', gene_list)
    logger.info('Searching on RCSB Server...')

    for gene in gene_list:
        gene = gene.strip()
        logger.info('\n[Searching %s]' % gene)

        query = attrs.rcsb_entity_source_organism.rcsb_gene_name.value == gene
        results = set(query())
        if not results:
            results = None

        logger.info('%s Search Results:' % gene, results)

        with open(gene + '.txt','w') as file:
            for id in results:
                file.write(id + '\n')
        
        logger.info('Search Complete.')