import csv
import json
import logging
import os
import re

logger = logging.getLogger('pyCADD.query.base')

import pandas as pd
from pyCADD.query import core
from pyCADD.utils.getinfo import get_project_dir
from pyCADD.utils.tool import mkdirs


class Query:
    '''
    与RCSB PDB服务器交互 获取基因相关晶体信息
    '''

    def __init__(self) -> None:

        self.gene_file = ''                     # 全部基因列表文件
        self.gene_list = []                     # 全部基因列表
        self.gene_dic = {}                      # 全部基因名称对照字典
        self.project_dir = get_project_dir()
        mkdirs(self.required_dirs)
    
    @property
    def query_dir(self):
        return self.project_dir + '/query'

    @property
    def pdblist_dir(self):
        return self.project_dir + '/pdblist'

    @property
    def json_dir(self):
        return self.project_dir + '/json'

    @property
    def info_dir(self):
        return self.project_dir + '/crystal_info'

    @property
    def affinity_dir(self):
        return self.project_dir + '/affinity'
    
    @property 
    def required_dirs(self):
        return [
            self.query_dir,
            self.pdblist_dir,
            self.json_dir,
            self.info_dir,
            self.affinity_dir
        ]

    def read_input(self):
        '''
        读取用户输入的需要查询的基因
        '''

        self.gene_dic = core.get_input()
        self.gene_list = list(self.gene_dic.keys())
        logger.info('Gene list: %s' % ''.join((gene + '(' + self.gene_dic[gene] + ') ') for gene in self.gene_list))

    def search(self, gene_list=[]):
        '''
        在RCSB服务器搜索列表中的基因相关晶体PDBID 存储为文本文件

        Parameters
        ----------
        gene_list : list
            要搜索的全部基因列表
        '''
        os.chdir(self.pdblist_dir)
        core.search_rcsb(self.gene_list)
        os.chdir(self.project_dir)

    def query(self):
        '''
        搜索晶体信息
        '''
        logger.info('Query Crystals Info from RCSB Server...')
        core.query_rcsb(self.gene_list)
        logger.info('Query Complete.')


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
                    if gene in line[x-1].upper().split(','):        # 判断哪一条链是目标基因的蛋白
                        match_chain = line[x]

                for x in ligchain_index:
                    if line[x-1] in ignore_list:
                        continue
                    if line[x] and line[x] in match_chain:  # 判断配体是否是结合目标基因蛋白的配体
                        inputfile.write('%s,%s,%s\n' % (pdb, line[x - 1], _flag))
                        

            f.close()
            inputfile.close()
