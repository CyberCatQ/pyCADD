import os
import sys
from numpy import index_exp
import openpyxl
import re
import pandas as pd
from openpyxl.styles import Font

root_path = os.path.abspath(os.path.dirname(__file__)).split('src')[0]  # 项目路径 绝对路径
result_path = root_path + 'lib/result/'                                 # 对接结果数据文件库
info_path = root_path + 'lib/info/'                                     # 晶体基本信息数据文件库
pdb_path = root_path.split('automatedMD')[0]                            # PDB项目绝对路径(如果有)
inputfile_path = root_path + 'lib/list/'                                # py4schrodinger输入文件库

font = Font(name='等线',size=12)

def result_merge(ligname=None):
    '''
    合并对接结果数据并输出到同一EXCEL中
    Parameter
    ----------
    ligname : str
        外源配体名称(如果有)
    '''
    withlig = ''
    if ligname:
        withlig = '_' + ligname
    mergefile = pdb_path + 'DOCK_FINAL_RESULTS%s.xlsx' % withlig

    resultfiles = os.popen("cd %s && ls | grep 'RESULTS%s.csv'" % (result_path, withlig)).read().splitlines()
    writer = pd.ExcelWriter(mergefile)
    for resultfile in resultfiles:
        gene = resultfile.split('_')[0]
        resultdata = pd.read_csv(result_path + resultfile)
        resultdata.to_excel(writer, gene, index=False)
    
    writer.save()

    wb = openpyxl.load_workbook(mergefile)
    for ws in wb:
        for row in ws.iter_rows():
            for cell in row:
                cell.font = font
    wb.save(mergefile)

def info_merge():
    '''
    合并晶体基础信息并分不同sheet输出到同一EXCEL中 并产生一张汇总sheet
    '''
    mergefile = pdb_path + 'CRYSTALS_INFO.xlsx'
    infofiles = os.popen('cd %s && ls' % info_path).read().splitlines()
    writer = pd.ExcelWriter(mergefile)
    data_list = []  # 汇总数据表

    for infofile in infofiles:
        gene = infofile.split('.')[0]
        data = pd.read_csv(info_path + infofile)
        data.to_excel(writer, gene, index=False)

        data['GENE'] = gene
        data_list.append(data)

    data_total = pd.concat(data_list)
    data_total.to_excel(writer, 'TOTAL', index=False)
    writer.save()

    wb = openpyxl.load_workbook(mergefile)
    for ws in wb:
        for row in ws.iter_rows():
            for cell in row:
                cell.font = font
    wb.save(mergefile)

def process_listfile():
    '''
    解析已存在的py4schrodinger输入文件 确保配体与目标基因的蛋白关联

    Return
    ---------
    dic
        { GENE: ['PDBID1,ligname1,agonist', 'PDBID2,ligname2,', ...] }
    '''
    list_dic = {}
    inputfiles = os.popen('cd %s && ls' % inputfile_path).read().splitlines()
    for inputfile in inputfiles:
        gene = inputfile.split('.')[0]
        with open(inputfile_path + inputfile) as f:
            cps = f.read().splitlines()
        list_dic[gene] = cps

    return list_dic

def conform_classify():
    '''
    晶体构象分类(Agonist/Antagonist)
    '''
    conform_dic = process_listfile()

    agonist_match = re.compile('agonist', re.IGNORECASE)
    antagonist_match = re.compile('antagonist', re.IGNORECASE)
    agonist_list = []
    antagonist_list = []
    agonist_ref = []
    antagonist_ref = []

    classifyfile = pdb_path + 'COMFORM_CLASSIFY.xlsx'
    infofiles = os.popen('cd %s && ls' % info_path).read().splitlines()
    writer = pd.ExcelWriter(classifyfile)

    for infofile in infofiles:
        gene = infofile.split('.')[0]
        data = pd.read_csv(info_path + infofile)
        agonist_ID = []
        antagonist_ID = []

        for cp in conform_dic[gene]:                # 根据文件中的flag识别为激动剂或拮抗剂
            _flag = cp.split(',')[2]
            _pdb = cp.split(',')[0]
            _lig = cp.split(',')[1]

            if re.search(antagonist_match, _flag):
                antagonist_ID.append(_pdb)
            elif re.search(agonist_match, _flag):
                agonist_ID.append(_pdb)

        # 解析CSV
        for index, row in data.iterrows():
            if row['PDBID'] in antagonist_ID:      # 拮抗剂命中
                antagonist_list.append({'Gene Name' : gene, 'PDB ID' : row['PDBID'],'Conformation' : 'Antagonist','Title' : row['title'], 'Reference' : row['reference'], 'DOI' : row['DOI']})
                antagonist_ref.append(row['DOI'])
            elif row['PDBID'] in agonist_ID:          # 激动剂命中
                agonist_list.append({'Gene Name' : gene, 'PDB ID' : row['PDBID'],'Conformation' : 'Agonist','Title' : row['title'], 'Reference' : row['reference'], 'DOI' : row['DOI']})
                agonist_ref.append(row['DOI'])

    agonist_df = pd.DataFrame(agonist_list)
    antagonist_df = pd.DataFrame(antagonist_list)
    agonist_df.to_excel(writer, 'Agonist', index=False)
    antagonist_df.to_excel(writer, 'Antagonist', index=False)

    # 文献去重
    antagonist_ref = set(antagonist_ref)
    agonist_ref = set(agonist_ref)
    antagonist_ref_df = pd.DataFrame(antagonist_ref)
    agonist_ref_df = pd.DataFrame(agonist_ref)
    agonist_ref_df.to_excel(writer, 'Agonist_Reference', index=False)
    agonist_ref_df.to_csv(pdb_path + 'Agonist_Ref.csv', index=False)
    antagonist_ref_df.to_excel(writer, 'Antagonist_Reference', index=False)
    antagonist_ref_df.to_csv(pdb_path + 'Antagonist_Ref.csv', index=False)
    writer.save()
    
    wb = openpyxl.load_workbook(classifyfile)
    for ws in wb:
        for row in ws.iter_rows():
            for cell in row:
                cell.font = font
    wb.save(classifyfile)

def main():
    print('''
    Script for Merging Data after Docking.

1. Merge the basic information of all crystals in automatedMD/lib/info
2. Merge all docking results data in automatedMD/lib/result
3. Merge all docking with exogenous ligand results data in automatedMD/lib/result
4. Classfify all crystals via title or reference

0. Exit
''')
    print('Enter the code for the operation to be performed:')
    _flag = input()
    if _flag == '1':
        info_merge()
    elif _flag == '2':
        result_merge()
    elif _flag == '3':
        ligname = input('Enter the ligand name:')
        result_merge(ligname)
    elif _flag == '4':
        conform_classify()
    else:
        sys.exit()
    
if __name__ == '__main__':
    main()