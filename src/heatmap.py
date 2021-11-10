import os
import re
from datetime import datetime

import matplotlib.pyplot as plt
import openpyxl
import pandas as pd
from pandas.core.series import Series
import seaborn as sns
from numpy import exp
from openpyxl.styles import Font

root_path = os.path.abspath(os.path.dirname(
    __file__)).split('src')[0]  # 项目路径 绝对路径
# PDB项目绝对路径(如果有)
pdb_path = root_path.split('automatedMD')[0]
final_result_path = pdb_path + 'final_results/'

now = datetime.now()
date = str(now.year).rjust(2, '0') + \
    str(now.month).rjust(2, '0') + str(now.day).rjust(2, '0')
font = Font(name='等线', size=12)
gene_dict = {'NR1A1': 'TR_alpha',
             'NR1A2': 'TR_beta',
             'NR1B1': 'RAR_alpha',
             'NR1B2': 'RAR_beta',
             'NR1B3': 'RAR_gamma',
             'NR1C1': 'PPAR_alpha',
             'NR1C2': 'PPAR_beta',
             'NR1C3': 'PPAR_gamma',
             'NR1D1': 'Rev-erb',
             'NR1F1': 'ROR_alpha',
             'NR1F2': 'ROR_beta',
             'NR1F3': 'ROR_gamma',
             'NR1H3': 'LXR_alpha',
             'NR1H2': 'LXR_beta',
             'NR1H4': 'FXR_alpha',
             'NR1H5': 'FXR_beta',
             'NR1I1': 'VDR',
             'NR1I2': 'PXR',
             'NR1I3': 'CAR',
             'NR2A1': 'HNF4_alpha',
             'NR2A2': 'HNF4_gamma',
             'NR2B1': 'RXR_alpha',
             'NR2B2': 'RXR_beta',
             'NR2B3': 'RXR_gamma',
             'NR2C1': 'TR2',
             'NR2C2': 'TR4',
             'NR2E2': 'TLL',
             'NR2E3': 'PNR',
             'NR2F1': 'COUP-TFI',
             'NR2F2': 'COUP-TFII',
             'NR2F6': 'EAR2',
             'NR3A1': 'ER_alpha',
             'NR3A2': 'ER_beta',
             'NR3B1': 'ERR_alpha',
             'NR3B2': 'ERR_beta',
             'NR3B3': 'ERR_gamma',
             'NR3C1': 'GR',
             'NR3C2': 'MR',
             'NR3C3': 'PR',
             'NR3C4': 'AR',
             'NR4A1': 'NGFIB',
             'NR4A2': 'NURR1',
             'NR4A3': 'NOR1',
             'NR5A1': 'SF1',
             'NR5A2': 'LRH1',
             'NR6A1': 'GCNF',
             'NR0B1': 'DAX1',
             'NR0B2': 'SHP'}
gene_list = gene_dict.keys()
gene_name = pd.Series(gene_dict, name='Abbreviation')
ignore_gene_list = ['NR5A1'] # ['NR1F1','NR2A1','NR1I3','NR3B1','NR1F2']
ignore_abbre_list = ['SF1'] # ['ROR_alpha','HNF4_alpha','CAR','ERR_alpha','ROR_beta']

def prefix_generate(date, module_name, precision) -> str:
    '''
    根据日期、模块和精度生成此次计算的结果和热力图文件名称前缀

    Parameters
    ----------
    date : str
        日期
    module_name : str
        生成数据的算法模块名
    precision : str
        结果来源的精度

    Return
    ----------
    str
        文件名前缀
    '''
    i = 1
    while True:
        file_prefix = date + '_' + \
            str(i).rjust(2, '0') + '_' + module_name + '_' + precision
        img_file_path = pdb_path + file_prefix + '.png'
        excel_file_path = pdb_path + file_prefix + '.xlsx'
        if os.path.exists(img_file_path) or os.path.exists(excel_file_path):
            i += 1
        else:
            break

    return file_prefix


def gen_heatmap(data, file_path) -> None:
    '''
    生成热力图并保存

    Parameters
    ----------
    data : pandas Dataframe
        生成热力图的DF数据对象
    file_path:
        热力图保存路径
    '''
    sns.set_context({'figure.figsize': (20, 20)})
    sns.heatmap(data=data, square=True, cmap="RdBu_r", annot=True,
                fmt=".0f", linewidths=0.1, vmin=0, vmax=100)
    plt.savefig(file_path, bbox_inches='tight')


class Dataprocessor:
    '''
    读取 处理和加工数据
    '''

    def __init__(self, file_path) -> None:
        file_path = os.path.abspath(file_path)
        data = pd.read_excel(file_path)
        self.origin_data = data[data['rmsd'].notnull()]      # 共结晶配体对接数据
        self.exligand_data = data[data['rmsd'].isnull()]     # 外源配体对接数据
        # 外源配体名称
        self.ex_ligand = re.search(
            r'(?<=FINAL_RESULTS_)[0-9a-zA-Z]+(?=_)', file_path).group()  
        # 各基因成员总数
        self.gene_count = self.origin_data['Gene Name'].value_counts()
        # 对接分数统计数据
        self.ds_data = pd.read_excel(file_path.split(
            '.')[0] + '_Pivot_table_Docking_Score.xlsx', sheet_name=self.ex_ligand)  
        # mmgbsa结合能统计数据
        self.mmgbsa_data = pd.read_excel(file_path.split('.')[
                                    0] + '_Pivot_table_MMGBSA_dG_Bind.xlsx', sheet_name=self.ex_ligand)  
        # Z-score
        self.z_score_ds = pd.read_excel(file_path.split('.')[
                                   0] + '_Pivot_table_Docking_Score.xlsx', sheet_name='Z-score', index_col='Abbreviation')['Z-score-combo']
        self.z_score_mmgbsa = pd.read_excel(file_path.split('.')[
                                       0] + '_Pivot_table_MMGBSA_dG_Bind.xlsx', sheet_name='Z-score', index_col='Abbreviation')['Z-score-combo']

    def top_count(self, by='Docking_Score') -> Series:
        '''
        对TOP100中含有各基因的成员计数 返回包含成员基因与计数值键值对的字典
        
        Parameters
        ----------
        by : str
            排序TOP时使用的数据列 默认为Docking_Score

        Return 
        ----------
        Series
            计数结果
        '''

        exligand_data = self.exligand_data
        exligand_data = exligand_data.sort_values(by=by, ascending=True).reset_index()
        # 所有计数初始化为0
        count = exligand_data[0:100]['Gene Name'].value_counts()            # 前100结果
        
        return count
    
    def calc_abs_score(self, by='Docking_Score') -> Series:
        '''
        计算所有基因由top_count()得到的计数占成员总数的比例值
        
        Return
        ----------
        Series
            比例值结果序列
        '''
        score = pd.Series(0.0, index = gene_list)
        gene_count = self.gene_count
        top_count = pd.Series(self.top_count(by))
        score = pd.Series(top_count / gene_count, name= by + '_ABS')

        return score
    
    def mmgbsa_penalty(self) -> Series:
        '''
        计算MMGBSA大于零的计数占成员总数的比例值

        Return
        ----------
        Series
            比例值结果序列
        '''
        gene_count = self.gene_count
        exligand_data = self.exligand_data
        penalty = exligand_data[exligand_data['MMGBSA_dG_Bind'] >= 0]['Gene Name'].value_counts()

        score = penalty / gene_count
        score = pd.Series(-score, name='Penalty')

        return score

    def calc_relative_score(self, by='Docking_Score') -> Series:
        '''
        计算考察变量的相对值 位于适当区间中的计数占成员总数的比例
        
        Parameter
        ----------
        by : str
            要考察的数据列

        Return
        ----------
        Series
            比例值结果序列
        '''

        if by == 'Docking_Score':
            _data = self.ds_data
        
        elif by == 'MMGBSA_dG_Bind':
            _data = self.mmgbsa_data

        else:
            raise RuntimeError('No property named %s' % by)

        gene_count = self.gene_count
        # 可变范围
        upper_limit = 1.3
        lower_limit = 0.7

        count = pd.Series(0, index=gene_list)
        for property in ['AVERAGE_All', 'Min_All', 'PDB_%s' % by]:
            tmp_df = _data[_data[property] >= lower_limit]
            tmp_df = tmp_df[tmp_df[property] <= upper_limit]
            count = count.add(tmp_df['Gene Name'].value_counts(), fill_value=0)
        
        score = pd.Series(count / (3 * gene_count), name=by+'_RELA')

        return score

    def calc_final_score(self):
        '''
        通识算法1 最终分数计算

        Return
        ----------
        tuple(Series, Dataframe)
            最终分数序列， 全部分数合并后的数据矩阵
        '''
        ligname = self.ex_ligand

        ds_abs_score = self.calc_abs_score(by='Docking_Score')
        mmgbsa_abs_score = self.calc_abs_score(by='MMGBSA_dG_Bind')
        mmgbsa_penalty = self.mmgbsa_penalty()
        ds_rela_score = self.calc_relative_score(by='Docking_Score')
        mmgbsa_rela_score = self.calc_relative_score(by='MMGBSA_dG_Bind')
        score_list = [ds_abs_score, mmgbsa_abs_score, mmgbsa_penalty, ds_rela_score, mmgbsa_rela_score]

        final_score = ds_abs_score
        for score in [mmgbsa_abs_score, mmgbsa_penalty, ds_rela_score, mmgbsa_rela_score]:
            final_score = final_score.add(score, fill_value=0)
        # final_score -= min(final_score)
        final_score = final_score / 4 * 100
        final_score.rename(ligname, inplace=True)

        score_list.append(final_score)

        df = pd.concat(score_list, axis=1)
        df.sort_index(inplace=True)
        for index, row in df.iterrows():
            df.loc[index, 'Abbreviation'] = gene_dict[index]
        df.set_index('Abbreviation', inplace=True)

        return final_score, df
    
    def calc_z_score(self):
        '''
        通识算法2 z-score计算

        Return
        ----------
        tuple(Series, Dataframe)
            Z-score 分数序列， 分数合并后的数据矩阵
        '''

        z_score_ds = self.z_score_ds.rename('Z_score_DS')
        z_score_mmgbsa = self.z_score_mmgbsa.rename('Z_score_MMGBSA')
        ligname = self.ex_ligand

        z_score_final = z_score_ds.add(z_score_mmgbsa, fill_value=0)
        z_score_final = 100 / exp(z_score_final)
        z_score_final.rename(ligname, inplace=True)

        df = pd.concat([z_score_ds, z_score_mmgbsa, z_score_final], axis=1)
        df.sort_index(inplace=True)

        return z_score_final, df


def main():
    '''
    程序入口
    '''
    _tmp = input('Enter all ligands need to be showed(split with comma):').split(',')
    precision = input('Enter the docking precision:').strip().upper()
    print('''
Calculation Method:

1. General Evaluation
2. Z-score
3. Machine Learning
    ''')
    flag = input('Enter the code of calculation method: ').strip()

    ligs = []
    for lig in _tmp:
        ligs.append(lig.strip())
    data = []

    if flag == '1':
        prefix = prefix_generate(date, 'GE', precision)
        excel_file_path = pdb_path + prefix + '.xlsx'
        writer = pd.ExcelWriter(excel_file_path)

        data = [gene_name]
        for lig in ligs:
            file_path = final_result_path + 'FINAL_RESULTS_%s_%s.xlsx' % (lig, precision)
            dataprocessor = Dataprocessor(file_path)
            result, df = dataprocessor.calc_final_score()
            # 移除异常的基因
            for gene in ignore_gene_list:
                result.drop(index=gene, inplace=True)
            for gene in ignore_abbre_list:
                df.drop(index=gene, inplace=True)

            data.append(result)
            df.to_excel(writer, sheet_name=lig)

    elif flag == '2':
        prefix = prefix_generate(date, 'Z', precision)
        excel_file_path = pdb_path + prefix + '.xlsx'
        writer = pd.ExcelWriter(excel_file_path)
        for lig in ligs:
            file_path = final_result_path + 'FINAL_RESULTS_%s_%s.xlsx' % (lig, precision)
            dataprocessor = Dataprocessor(file_path)
            result, df = dataprocessor.calc_z_score()
            # 移除异常的基因
            for gene in ignore_abbre_list:
                result.drop(index=gene, inplace=True)
                df.drop(index=gene, inplace=True)

            data.append(result)
            df.to_excel(writer, sheet_name=lig)
    else:
        raise RuntimeError('Wrong Code, Exit.')
    
    data = pd.concat(data, axis=1)
    try:
        data.set_index('Abbreviation', inplace=True)
    except KeyError:
        pass
    data.sort_index(inplace=True)
    data.dropna(axis=0, how='all', inplace=True)
    data.to_excel(writer, sheet_name='TOTAL')
    
    writer.save()

    wb = openpyxl.load_workbook(excel_file_path)
    for ws in wb:
        for row in ws.iter_rows():
            for cell in row:
                cell.font = font
    wb.save(excel_file_path)

    heatmap_file_path = pdb_path + prefix + '.png'
    gen_heatmap(data, heatmap_file_path)

if __name__ == '__main__':
    main()




