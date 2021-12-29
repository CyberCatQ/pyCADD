from pyCADD.utils.check import check_file
from urllib import parse
import os
import logging
import requests
import json
from time import sleep

from pyCADD.utils.getinfo import get_project_dir
logger = logging.getLogger('pyCADD.query.core')
cwd = get_project_dir()

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

def _search_gene(gene:str):
    '''
    搜索基因相关PDB
    Parameter
    ---------
    gene : str
        查询基因名
    Return
    ---------
    list
        PDB ID 列表
    '''
    url = 'https://search.rcsb.org/rcsbsearch/v1/query'
    query = {
    "query":{
        "type" : "terminal",
        "service" : "text",
        "parameters":{
            'attribute': "rcsb_entity_source_organism.rcsb_gene_name.value",
            "operator" : "exact_match",
            "value" : gene
        }
    },
    "return_type": "entry"
            }
    data = json.dumps(query)

    response = requests.post(url, data)
    logger.debug('Searching %s status code: %s' % (gene, response.status_code))
    # 空结果情况
    if response.status_code != 200:
        return None

    results_text = json.loads(response.text)
    results = []

    for result_item in results_text['result_set']:
        results.append(result_item['identifier'])
    
    return results
    
def search_rcsb(gene_list:list):
    '''
    在RCSB服务器搜索列表中的基因相关晶体PDBID 存储为文本文件

    Parameters
    ----------
    gene_list : list
        要搜索的全部基因列表
    '''


    logger.debug('\nAll Genes: %s' % gene_list)
    logger.info('Searching on RCSB Server...')
    os.chdir(cwd + '/pdblist')

    for gene in gene_list:
        gene = gene.strip()
        logger.info('[Searching %s]' % gene)
        results = _search_gene(gene)

        if results:
            with open(gene + '.txt','w') as file:
                for id in results:
                    file.write(id + '\n')

        logger.info('%s Search Results: %s' % (gene, results))
    
    logger.info('Search Complete.')
    os.chdir(cwd)

def __save_data(response, gene):
    '''
    发送请求内容 存储从服务器返回的信息为json文件

    Parameters
    ----------
    response : requsets.Response
        响应数据
    gene : str
        用以命名json文件的基因名称
    '''
        
    with open(gene + '.json', 'w') as f:                    # 存储为json
        f.write(response.text)

def query_rcsb(gene_list=[]):
    '''
    向RCSB服务器提交请求 获取PDB晶体信息
    请求晶体信息内容:
        PDB ID
        晶体标题
        晶体相关文献
        文献作者
        晶体配体亲和力实验数据(如果有)
        配体总数量
        聚体(肽链)总数量
        配体名称
        配体RCSB ID
        配体结合肽链编号
        配体SMILES式
        肽链名称
        肽链源基因名
        肽链编号
        
    Parameters
    ----------
    gene_list: list
        要请求的基因列表
    '''

    pdbids = ''
    url = 'https://data.rcsb.org/graphql?query='                            # API接口

    os.chdir(cwd + '/json')

    for gene in gene_list:
        pdbids = ''
        gene = gene.strip()
        
        if check_file(gene + '.json'):
            logger.debug('%s Crystals Information Already Exists. Skip.' % gene)
            continue
        logger.info('Getting Crystals Information: %s' % gene)

        with open(cwd + '/pdblist/%s.txt' % gene, 'r') as f:
            pdbs = f.read().splitlines()
            if pdbs:                                                        # 将所有PDBID 以逗号分割并相连 组成请求字符串
                for i in pdbs:
                    pdbids = pdbids + "\"" + i + "\"" + ','
                pdbids = pdbids.strip(',')

                # 请求内容 与解析部分关联 谨慎修改
                query = '''\
{
	entries(entry_ids: [%s])
	{
		rcsb_id
		rcsb_primary_citation {
			pdbx_database_id_DOI
			rcsb_authors
			title
		}
		rcsb_binding_affinity {
		comp_id
			value
			provenance_code
			type
			unit
		}
		rcsb_entry_container_identifiers {
			entry_id
		}
		rcsb_entry_info {
			deposited_nonpolymer_entity_instance_count
			deposited_polymer_entity_instance_count
		}
		struct {
			title
		}
		polymer_entities {
			entity_poly {
				rcsb_entity_polymer_type
			}
			rcsb_entity_source_organism {
				rcsb_gene_name {
					value
				}
			}
			rcsb_polymer_entity {
				pdbx_description
			}
			rcsb_polymer_entity_container_identifiers {
				auth_asym_ids
			}
		}
		nonpolymer_entities {
			nonpolymer_comp {
				chem_comp {
					id
					name
				}
				rcsb_chem_comp_descriptor {
					SMILES
				}
			}
			rcsb_nonpolymer_entity_container_identifiers {
				auth_asym_ids
			}
		}
	}
}
    ''' % pdbids
                query = parse.quote(query)              # 将请求内容以url编码
                request = url + query                   # 组装请求url
                response = requests.get(request)
                __save_data(response, gene)         # 发送请求 保存返回数据为json文件
                sleep(2)
            else:
                logger.warning('%s.txt has no crystal.' % gene)

    os.chdir(cwd)
