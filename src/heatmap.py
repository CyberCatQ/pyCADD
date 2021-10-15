import pandas as pd
import os
import re
import seaborn as sns
from datetime import datetime
import matplotlib.pyplot as plt

root_path = os.path.abspath(os.path.dirname(__file__)).split('src')[0]  # 项目路径 绝对路径
pdb_path = root_path.split('automatedMD')[0]                            # PDB项目绝对路径(如果有)
now = datetime.now()
date = str(now.year) + str(now.month) + str(now.day)

writer = pd.ExcelWriter(pdb_path + '%s_Scores.xlsx' % date)
gene_dict = {'NR1A1': 'TR_alpha',
             'NR1A2': 'TR_beta',
             'NR1B1': 'RAR_alpha',
             'NR1B2': 'RAR_beta',
             'NR1B3': 'RAR_gamma',
             'NR1C1': 'PPAR_alpha',
             'NR1C2': 'PPAR_beta',
             'NR1C3': 'PPAR_gamma',
             'NR1F1': 'ROR_alpha',
             'NR1F2': 'ROR_beta',
             'NR1F3': 'ROR_gamma',
             'NR1H3': 'LXR_alpha',
             'NR1H2': 'LXR_beta',
             'NR1H4': 'FXR_alpha',
             'NR1I1': 'VDR',
             'NR1I2': 'PXR',
             'NR1I3': 'CAR',
             'NR2A1': 'HNF4_alpha',
             'NR2A2': 'HNF4_gamma',
             'NR2B1': 'RXR_alpha',
             'NR2B2': 'RXR_beta',
             'NR2B3': 'RXR_gamma',
             'NR3A1': 'ER_alpha',
             'NR3A2': 'ER_beta',
             'NR3B1': 'ERR_alpha',
             'NR3B3': 'ERR_gamma',
             'NR3C1': 'GR',
             'NR3C2': 'MR',
             'NR3C3': 'PR',
             'NR3C4': 'AR',
             'NR4A1': 'NGFIB',
             'NR4A2': 'NURR1',
             'NR5A1': 'SF1',
             'NR5A2': 'LRH1'}
gene_name = pd.Series(gene_dict, name='Abbreviation')

class Heatmap:
    '''
    识别情况热力图
    '''

    def __init__(self, file_path) -> None:

        file_path = os.path.abspath(file_path)
        raw_data = pd.read_excel(file_path)
        data = raw_data[raw_data['Docking_Score'].notnull()]
        origin_data = data[data['rmsd'].notnull()]
        exligand_data = data[data['rmsd'].isnull()]
        gene_count = origin_data['Gene Name'].value_counts().to_dict()
        ex_ligand = re.search(r'(?<=FINAL_RESULTS_)[0-9A-Z]+(?=_)', file_path).group()
        ds_data = pd.read_excel(file_path.split('.')[0] + '_Pivot_table_Docking_Score.xlsx', sheet_name=ex_ligand)
        mmgbsa_data = pd.read_excel(file_path.split('.')[0] + '_Pivot_table_MMGBSA_dG_Bind.xlsx', sheet_name=ex_ligand)

        self.raw_data = raw_data
        self.data = data
        self.origin_data = origin_data
        self.exligand_data = exligand_data
        self.gene_count = gene_count
        self.gene_list = list(gene_count.keys())
        self.ds_data = ds_data
        self.mmgbsa_data = mmgbsa_data
        self.ligname = file_path.split('/')[-1].split('_')[2]

    def _top_count(self, by='Docking_Score'):

        exligand_data = self.exligand_data
        exligand_data = exligand_data.sort_values(by=by, ascending=True).reset_index()
        top_count = exligand_data[0:100]['Gene Name'].value_counts().to_dict()              # 前100结果
        
        return top_count

    def calc_abs_score(self, dic):
        
        score_dict = {}
        gene_count = self.gene_count
        gene_list = self.gene_list

        for gene in gene_list:
            score_dict[gene] = 0
            try:
                score_dict[gene] = dic[gene] / gene_count[gene]
            except KeyError:
                continue

        return score_dict
    
    def hit_percent(self):

        exligand_data = self.exligand_data
        gene_count = self.gene_count
        gene_list = self.gene_list
        hit_dict = {}
        
        hit = exligand_data['Gene Name'].value_counts()                                     # 各基因的命中数
        for gene in gene_list:
            hit_dict[gene] = 0
            try:
                hit_dict[gene] = hit[gene] / gene_count[gene]
            except KeyError:
                continue
            
        return hit_dict

    def mmgbsa_punish(self):
        
        gene_count = self.gene_count
        exligand_data = self.exligand_data
        punish_count = {}
        punish_dict = {}

        for gene in gene_count.keys():
            punish_count[gene] = 0

        for index, row in exligand_data.iterrows():
            if row['MMGBSA_dG_Bind'] > 0:
                punish_count[row['Gene Name']] += 1                                         # MMGBSA大于0的次数
        
        for gene in gene_count.keys():
            punish_dict[gene] = punish_count[gene] / gene_count[gene]                       # 大于0的数目占该基因总成员百分比

        return punish_dict
    
    def calc_relative_score(self, by='Docking_Score'):

        if by == 'Docking_Score':
            _data = self.ds_data
        
        elif by == 'MMGBSA_dG_Bind':
            _data = self.mmgbsa_data

        else:
            raise RuntimeError('No property named %s' % by)

        gene_count = self.gene_count
        rela_calc = {}
        rela_score_dict = {}

        for gene in gene_count.keys():
            rela_calc[gene] = 0
            rela_score_dict[gene] = 0

        for index, row in _data.iterrows():
            if row['AVERAGE_All'] >= 0.8 and row['AVERAGE_All'] <= 1.2:
                rela_calc[row['Gene Name']] += 1
            if row['Min_All'] >= 0.8 and row['Min_All'] <= 1.2:
                rela_calc[row['Gene Name']] += 1
            if row['PDB_%s' % by] >= 0.8 and row['PDB_%s' % by] <= 1.2:
                rela_calc[row['Gene Name']] += 1

        for gene in rela_calc.keys():
            rela_score_dict[gene] = rela_calc[gene] / (3 * gene_count[gene])

        return rela_score_dict
    
    def calc_final_scroe(self):
        
        ligname = self.ligname
        exligand_data = self.exligand_data
        gene_list = exligand_data['Gene Name'].value_counts().keys()

        ds_abs_score = self.calc_abs_score(self._top_count(by='Docking_Score'))
        # print('DS ABS SCORE:', ds_abs_score)
        mmgbsa_abs_score = self.calc_abs_score(self._top_count(by='MMGBSA_dG_Bind'))
        # print('\nMMGBSA ABS SCORE:', mmgbsa_abs_score)

        hit_percent = self.hit_percent()
        # print('\nHit Percent:', hit_percent)
        mmgbsa_punish = self.mmgbsa_punish()
        # print('\nMMGBSA Punish:', mmgbsa_punish)

        ds_rela_score = self.calc_relative_score(by='Docking_Score')
        # print('\nDS RELATIVE SCORE:', ds_rela_score)
        mmgbsa_rela_score = self.calc_relative_score(by='MMGBSA_dG_Bind')
        # print('\nMMGBSA RELATIVE SCORE:', mmgbsa_rela_score)

        final_score = {}
        for gene in gene_list:
            final_score[gene] = (ds_abs_score[gene] + mmgbsa_abs_score[gene] + hit_percent[gene] - mmgbsa_punish[gene] + ds_rela_score[gene] + mmgbsa_rela_score[gene]) / 5 * 100
        # print('\nFinal Score:', final_score)
        
        series1 = pd.Series(ds_abs_score, name='DS_ABS')
        series2 = pd.Series(mmgbsa_abs_score, name='MMGBSA_ABS')
        series3 = pd.Series(hit_percent, name='HIT')
        series4 = pd.Series(mmgbsa_punish, name='PUNISH')
        series5 = pd.Series(ds_rela_score, name='DS_RELA')
        series6 = pd.Series(mmgbsa_rela_score, name='MMGBSA_RELA')
        series7 = pd.Series(final_score, name=self.ligname)

        df = pd.concat([series1,series2,series3,series4,series5,series6,series7], axis=1)
        
        df.sort_index(inplace=True)
        df.to_excel(writer, sheet_name=ligname)

        return series7

    @staticmethod
    def gen_heatmap(df):
        
        sns.set_context({'figure.figsize':(16,16)})
        sns.heatmap(data=df, square=True,cmap="RdBu_r", annot=True, fmt=".0f" ,linewidths=0.1)
        plt.show()
        plt.savefig(pdb_path + '%s_heatmap.png' % date, bbox_inches='tight')

def main():
    ligs = []
    _tmp = input('Enter all ligands need to be showed(split with comma):').split(',')
    precision = input('Enter the docking precision:').strip().upper()
    for lig in _tmp:
        ligs.append(lig.strip())
    
    data = [gene_name]
    for lig in ligs:
        file_path = pdb_path + 'FINAL_RESULTS_%s_%s.xlsx' % (lig, precision)
        heatmap = Heatmap(file_path)
        result = heatmap.calc_final_scroe()
        data.append(result)

    data = pd.concat(data, axis=1)
    data.sort_index(inplace=True)
    data.set_index('Abbreviation', inplace=True)
    data.to_excel(writer, sheet_name='TOTAL')
    writer.save()
    heatmap.gen_heatmap(data)
    
if __name__ == '__main__':
    main()