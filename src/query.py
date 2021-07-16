import csv
import json
import os
import sys
import re
from time import sleep
from urllib import parse

import pandas as pd
import requests
from rcsbsearch import rcsb_attributes as attrs

root_path = os.path.abspath(os.path.dirname(__file__)).split('src')[0]  # 根路径 绝对路径
lib_path = root_path + 'lib/'                                           # 库文件夹路径
doc_path = root_path + 'doc/'                                           # 文档文件夹路径

pdblist_dir_path = lib_path + 'pdbid/'                              # 记录基因PDBID的列表文件库
inputfile_dir_path = lib_path + 'list/'                             # py4schrodinger输入文件库
json_dir_path = lib_path + 'json/'                                  # 服务器返回的基因原始信息json库
info_dir_path = lib_path + 'info/'                                  # 基因晶体信息转换后的csv格式文件库
bdaffinity_dir_path = lib_path + 'BDaffinity/'                      # 亲和力实验数据库

def make_dir(path):
    try:
        os.makedirs(path)
    except FileExistsError:
        pass

make_dir(pdblist_dir_path)
make_dir(inputfile_dir_path)
make_dir(json_dir_path)
make_dir(info_dir_path)
make_dir(bdaffinity_dir_path)

class Query:
    '''
与RCSB PDB服务器交互 获取基因相关晶体信息
    '''

    def __init__(self) -> None:

        self.gene_file = ''         # 全部基因列表文件
        self.gene_list = []         # 全部基因列表
        self.gene_dic = {}          # 全部基因名称对照字典

    def _precess_inputfile(self):
        '''
        解析输入文件
        '''
        gene_file = sys.argv[1]      # 从命令行解析提供的基因列表文件

        try:
            dic = {}
            with open(gene_file,'r') as f:
                gene_infos = f.readlines()
                for gene_info in gene_infos:
                    gene_info = gene_info.strip()
                    gene = gene_info.split(',')[0]
                    name = gene_info.split(',')[1]
                    dic[gene] = name
                self.gene_dic = dic
                self.gene_list = list(dic.keys())
                self.gene_file = os.path.abspath(gene_file)
            print('\nGene List File: %s ' % gene_file)
        except FileNotFoundError:
            print('Error: File %s not Found' % gene_file)
            raise

    def search(self, gene_list=[]):
        '''
        在RCSB服务器搜索列表中的基因相关晶体PDBID 存储为文本文件

        Parameters
        ----------
        gene_list : list
            要搜索的全部基因列表
        '''
        os.chdir(pdblist_dir_path)
        if not gene_list:
            gene_list = self.gene_list

        print('\nAll Genes:', gene_list)
        print('Searching on RCSB Server...')

        for gene in gene_list:
            gene = gene.strip()
            print('\n[Searching %s]' % gene)
            query = attrs.rcsb_entity_source_organism.rcsb_gene_name.value == gene
            results = set(query())
            if results:
                print('%s Search Results:' % gene, results)
            else:
                print('%s Search Results:' % gene, None)
            with open(gene + '.txt','w') as file:
                for id in results:
                    file.write(id + '\n')
        
        print('Search Complete.')

    def __save_data(self, request, gene):
        '''
        发送请求内容 存储从服务器返回的信息为json文件

        Parameters
        ----------
        request : str
            请求内容
        gene : str
            用以命名json文件的基因名称
        '''
        
        response = requests.get(request)                        # request get方法发送请求 获取response内容
        json_data = response.text                               # response内容以json解析
        with open(json_dir_path + gene + '.json','w') as f:     # 存储为json
            f.write(json_data)

    def query(self, gene_list=[]):
        '''
        向RCSB服务器提交请求 获取PDB晶体信息
        请求晶体信息内容：
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
        print('Query Crystals Info from RCSB Server...')
        if not gene_list:
            gene_list = self.gene_list                                          # 读取所有列表文件名称

        pdbids = ''
        url = 'https://data.rcsb.org/graphql?query='                            # API接口
        for gene in gene_list:
            pdbids = ''
            gene = gene.strip()
            
            if os.path.exists(json_dir_path + gene + '.json'):
                print('%s Crystals Information Already Exists. Skip.' % gene)
                continue
            print('[Get %s Crystals Information]' % gene)

            with open(pdblist_dir_path + gene + '.txt', 'r') as f:
                pdbs = f.readlines()
                if pdbs:                                                        # 将所有PDBID 以逗号分割并相连 组成请求字符串
                    for i in pdbs:
                        pdbids = pdbids + "\"" + i.strip() + "\"" + ','
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

                    self.__save_data(request, gene)         # 发送请求 保存返回数据为json文件
                    sleep(2)
                else:
                    print('WARNING: %s.txt has no crystal.' % gene)
        print('\nQuery Complete.')


    def convert_csv(self, gene_list=[]):
        '''
        转换json文件为csv
        '''

        '''
        if not gene_dic:
            gene_dic = self.gene_dic
        '''
        if not gene_list:
            gene_list = self.gene_list

        for gene in gene_list:
            _temp_value = []
            gene = gene.strip()

            try:
                jsfile = open(json_dir_path + gene + '.json', 'r')
            except FileNotFoundError:
                continue
            js_data = json.load(jsfile)
            lis_data = js_data['data']['entries']
            polymer_num = 0
            nonpolymer_num = 0

            for dic_pdb in lis_data:

                if polymer_num < len(dic_pdb['polymer_entities']):
                    polymer_num = len(dic_pdb['polymer_entities'])

                d = {}
                d['PDBID'] = dic_pdb['rcsb_id']  # 字符串
                d['nonpolymer_entity'] = dic_pdb['rcsb_entry_info']['deposited_nonpolymer_entity_instance_count']  # 字符串
                d['polymer_entity_instance'] = dic_pdb['rcsb_entry_info']['deposited_polymer_entity_instance_count']  # 字符串
                d['title'] = dic_pdb['struct']['title']  # 字符串
                d['reference'] = dic_pdb['rcsb_primary_citation']['title']  # 字符串
                d['authors'] = ''.join(
                    i + ' ' for i in dic_pdb['rcsb_primary_citation']['rcsb_authors'])  # 字符串
                d['DOI'] = dic_pdb['rcsb_primary_citation']['pdbx_database_id_DOI']  # 字符串
                d['polymer_entities'] = dic_pdb['polymer_entities']  # 多维列表字典
                d['nonpolymer_entities'] = dic_pdb['nonpolymer_entities']  # 多维列表字典
                d['binding_affinity'] = dic_pdb['rcsb_binding_affinity']  # 列表字典

                j = 0

                for i in d['polymer_entities']:

                    d['polymer_entities_' +
                        str(j)] = d['polymer_entities'][j]  # 多聚体对象
                    d['polymer_entities_' + str(j) + '_type'] = d['polymer_entities_' + str(
                        j)]['entity_poly']['rcsb_entity_polymer_type']  # 多聚体类型
                    d['polymer_entities_' + str(j) + '_name'] = d['polymer_entities_' + str(
                        j)]['rcsb_polymer_entity']['pdbx_description']  # 多聚体名称

                    if d['polymer_entities_' + str(j)]['rcsb_entity_source_organism']:
                        # 多聚体源基因嵌套字典列表
                        temp_list = d['polymer_entities_' +
                                    str(j)]['rcsb_entity_source_organism'][0]['rcsb_gene_name']
                        if temp_list:
                            s = ''
                            for di in temp_list:  # 每一个di是一个字典
                                s = s + di['value'] + ','
                        # 提取并拼接基因名称 用逗号相连存为字符串
                        d['polymer_entities_' + str(j) + '_gene'] = s.strip(',')

                    temp_list = d['polymer_entities_' +
                                str(j)]['rcsb_polymer_entity_container_identifiers']['auth_asym_ids']  # 多聚体所在链列表
                    s = ''
                    for t in temp_list:
                        s = s + t + ','
                    # 提取并拼接链列表 用逗号相连存为字符串
                    d['polymer_entities_' + str(j) + '_chain'] = s.strip(',')
                    del d['polymer_entities_' + str(j)]

                    j += 1
                del d['polymer_entities']

                if d['nonpolymer_entities']:
                    if nonpolymer_num < len(dic_pdb['nonpolymer_entities']):
                        nonpolymer_num = len(dic_pdb['nonpolymer_entities'])
                    j = 0

                    for i in d['nonpolymer_entities']:
                        s = ''
                        d['nonpolymer_entities_' +
                        str(j)] = d['nonpolymer_entities'][j]  # 小分子对象
                        d['nonpolymer_entities_' + str(j) + '_ligname'] = d['nonpolymer_entities_' + str(
                            j)]['nonpolymer_comp']['chem_comp']['name']  # 小分子名称
                        d['nonpolymer_entities_' + str(j) + '_ligid'] = d['nonpolymer_entities_' + str(
                            j)]['nonpolymer_comp']['chem_comp']['id']  # 小分子ID
                        d['nonpolymer_entities_' + str(j) + '_SMILES'] = d['nonpolymer_entities_' + str(
                            j)]['nonpolymer_comp']['rcsb_chem_comp_descriptor']['SMILES']

                        temp_list = d['nonpolymer_entities_' + str(
                            j)]['rcsb_nonpolymer_entity_container_identifiers']['auth_asym_ids']  # 小分子所在链列表
                        for t in temp_list:
                            s = s + t + ','
                        # 提取并拼接链列表 用逗号相连存为字符串
                        d['nonpolymer_entities_' +
                            str(j) + '_ligchain'] = s.strip(',')

                        del d['nonpolymer_entities_' + str(j)]
                        j += 1
                del d['nonpolymer_entities']

                if d['binding_affinity']:
                    with open(lib_path + 'BDaffinity/' + d['PDBID'] + '_binding_affinity.csv', 'w') as f:
                        df = pd.DataFrame(d['binding_affinity'], columns=[
                                        'comp_id', 'value', 'provenance_code', 'type', 'unit'])
                        df.sort_values(
                            by=['comp_id', 'provenance_code', 'type'], ascending=True, inplace=True)
                        df = df.reset_index(drop=True)
                        df.to_csv(f)

                del d['binding_affinity']

                _temp_value.append(d)

            keys = ['PDBID', 'title', 'polymer_entity_instance', 'nonpolymer_entity']
            for a in range(polymer_num):
                keys.append('polymer_entities_%s_name' % str(a))
                keys.append('polymer_entities_%s_type' % str(a))
                keys.append('polymer_entities_%s_gene' % str(a))
                keys.append('polymer_entities_%s_chain' % str(a))
            for b in range(nonpolymer_num):
                keys.append('nonpolymer_entities_%s_ligname' % str(b))
                keys.append('nonpolymer_entities_%s_ligid' % str(b))
                keys.append('nonpolymer_entities_%s_ligchain' % str(b))
                keys.append('nonpolymer_entities_%s_SMILES' % str(b))
                # keys.append('nonpolymer_entities_%s_activity' % str(b))

            keys.append('reference')
            keys.append('DOI')
            keys.append('authors')

            with open(info_dir_path + gene + '.csv', 'w') as f:
                writer = csv.DictWriter(f, fieldnames=keys)
                writer.writeheader()
                writer.writerows(_temp_value)

    def make_inputfile(self, gene_list=[]):
        '''
        依据基因名称 自动识别晶体中对应链及附属配体 以此生成py4schrodinger.py可读取的输入文件
        '''
        def __try_index(lis, obj):
            try:
                index = lis.index(obj)
            except ValueError:
                return None
            else:
                return index

        ignore_list = ['EDO', 'DMS', 'IPA', 'TBY', 'ARS', 'EU', 'MG', 'IOD', 'ACT', 'CA', 'CAC', 'K',
               'PO4', 'BR', 'NO3', 'BCT', 'ZN', 'SO4', 'CL', 'NA', 'AU', 'GOL', 'NI', 'YT3']  # 非配体例外列表

        '''
        EDO: 1,2-ETHANEDIOL
        DMS: DIMETHYL SULFOXIDE
        IPA: ISOPROPYL ALCOHOL
        ARS: ARSENIC
        EU:  EUROPIUM ION
        MG:  MAGNESIUM ION
        IOD: IODIDE ION
        ACT: ACETATE ION
        CA:  CALCIUM ION
        CAC: CACODYLATE ION
        K:   POTASSIUM ION
        PO4: PHOSPHATE ION
        BR:  BROMIDE ION
        NO3: NITRATE ION
        BCT: BICARBONATE ION
        ZN:  ZINC ION
        SO4: SULFATE ION
        CL:  CHLORIDE ION
        NA:  SODIUM ION
        AU:  GOLD ION
        GOL: GLYCEROL
        NI:  NICKEL (II) ION
        YT3: YTTRIUM (III) ION
        '''

        agonist_match = re.compile('agonist', re.IGNORECASE)
        antagonist_match = re.compile('antagonist', re.IGNORECASE)

        if not gene_list:
            gene_list = self.gene_list

        for gene in gene_list:
            
            csv_file = info_dir_path + gene + '.csv'
            try:
                f = open(csv_file, 'r')
            except FileNotFoundError:
                continue
            inputfile = open(inputfile_dir_path + gene + '.txt', 'w')
            reader = csv.reader(f)

            line1 = next(reader)
            chain_index = []
            ligchain_index = []

            for n in range(5):
                index = __try_index(line1, 'polymer_entities_%s_chain' % str(n))
                if index:
                    chain_index.append(index)

            for n in range(8):
                index = __try_index(
                    line1, 'nonpolymer_entities_%s_ligchain' % str(n))
                if index:
                    ligchain_index.append(index)

            for line in reader:
                _flag = ''
                # 判断构象类型
                if re.search(antagonist_match, str(line)):
                    _flag = 'antagonist'
                elif re.search(agonist_match, str(line)):
                    _flag = 'agonist'

                pdb = line[0]
                for x in chain_index:
                    if gene in line[x-1].split(','):        # 判断哪一条链是目标基因的蛋白
                        match_chain = line[x]

                for x in ligchain_index:
                    if line[x-1] in ignore_list:
                        continue
                    if line[x] and line[x] in match_chain:  # 判断配体是否是结合目标基因蛋白的配体
                        inputfile.write('%s,%s,%s\n' % (pdb, line[x - 1], _flag))
                        

            f.close()
            inputfile.close()

def ui():
    '''
    用户界面
    '''

    print(''.center(80, '*'))
    print('Gene Query Module'.center(80))
    print(''.center(80, '*'))
    query = Query()
    query._precess_inputfile()

    while True:
        flag = input('''
Input Query Code:

1. Search Crystals' PDB ID of Gene on RCSB Server
2. Get Basic Info of Crystals with Query
3. Make Input Files for py4schrodinger.py

0. Exit

Code:''')

        if flag == '1':
            query.search()
        elif flag == '2':
            query.query()
            query.convert_csv()
        elif flag == '3':
            query.make_inputfile()

        elif flag == '0':
            break

def _usage():
    '''
Search Genes ralated crystals PDBID and query basic info from RCSB server.

Usage: query <genes_list_file>

Description
    genes_list_file           Text file, including all genes needed to be search and query

Example for geneslist file:
    NR1A1,TR_alpha
    NR1A2,TR_beta
    NR1B1,RAR_alpha
    NR1B2,RAR_beta
    NR1B3,RAR_gamma
    NR1C1,PPAR_alpha
    NR1C2,PPAR_beta
    NR1C3,PPAR_gamma
    '''

if __name__ == '__main__':
    if len(sys.argv) == 1:
        print(_usage.__doc__)
    elif sys.argv[1] == '-h' or sys.argv[1] == '--help':
            print(_usage.__doc__)
    else:
        ui()
        
        
    



