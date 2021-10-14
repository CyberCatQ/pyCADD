import pandas as pd
import os
import re

class Heatmap:
    '''
    识别情况热力图
    '''

    def __init__(self, file_path) -> None:

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

    def _top_count(self, by='Docking_Score'):

        exligand_data = self.exligand_data
        exligand_data = exligand_data.sort_values(by=by, ascending=True).reset_index()
        top_count = exligand_data[0:100]['Gene Name'].value_counts().to_dict()              # 前100结果
        
        return top_count

    def calc_abs_score(self, dic):
        
        score_dict = {}
        gene_list = self.gene_list

        for gene in gene_list:
            score_dict[gene] = 0
            try:
                score_dict[gene] = dic[gene] / sum(dic.values())
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
        ral_calc = {}
        ral_score_dict = {}

        for gene in gene_count.keys():
            ral_calc[gene] = 0
            ral_score_dict[gene] = 0

        for index, row in _data.iterrows():
            if row['AVERAGE_All'] >= 0.8 and row['AVERAGE_All'] <= 1.2:
                ral_calc[row['Gene Name']] += 1

        for gene in ral_calc.keys():
            ral_score_dict[gene] = ral_calc[gene] / gene_count[gene]

        return ral_score_dict
    
    def gen_heatmap(self):
        
        exligand_data = self.exligand_data
        gene_list = exligand_data['Gene Name'].value_counts().keys()

        ds_abs_score = self.calc_abs_score(self._top_count(by='Docking_Score'))
        print('DS ABS SCORE:', ds_abs_score)
        mmgbsa_abs_score = self.calc_abs_score(self._top_count(by='MMGBSA_dG_Bind'))
        print('MMGBSA ABS SCORE:', mmgbsa_abs_score)

        hit_percent = self.hit_percent()
        mmgbsa_punish = self.mmgbsa_punish()

        ds_ral_score = self.calc_relative_score(by='Docking_Score')
        print('DS RALATIVE SCORE:', ds_ral_score)
        mmgbsa_ral_score = self.calc_relative_score(by='MMGBSA_dG_Bind')
        print('MMGBSA RALATIVE SCORE:', mmgbsa_ral_score)

        final_score = {}
        for gene in gene_list:
            final_score[gene] = ds_abs_score[gene] + mmgbsa_abs_score[gene] + hit_percent[gene] - mmgbsa_punish[gene] + ds_ral_score[gene] + mmgbsa_ral_score[gene]

        print('Final Score:', final_score)

def main():
    heatmap = Heatmap(r'C:\Users\54415\Desktop\HeatMap\FINAL_RESULTS_BRL_SP.xlsx')
    heatmap.gen_heatmap()

if __name__ == '__main__':
    main()