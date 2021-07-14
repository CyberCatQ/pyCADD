import os
import sys
import openpyxl
import pandas as pd
from openpyxl.styles import Font

root_path = os.path.abspath(os.path.dirname(__file__)).split('src')[0]  # 项目路径 绝对路径
result_path = root_path + 'lib/result/'                                 # 对接结果数据文件库
info_path = root_path + 'lib/info/'                                      # 晶体基本信息数据文件库
pdb_path = root_path.split('automatedMD')[0]                            # PDB项目绝对路径(如果有)
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
    合并晶体基础信息并输出到同一EXCEL中
    '''
    mergefile = pdb_path + 'CRYSTALS_INFO.xlsx'
    infofiles = os.popen('cd %s && ls' % info_path).read().splitlines()
    writer = pd.ExcelWriter(mergefile)

    for infofile in infofiles:
        gene = infofile.split('.')[0]
        data = pd.read_csv(info_path + infofile)
        data.to_excel(writer, gene, index=False)
    writer.save()

    wb = openpyxl.load_workbook(mergefile)
    for ws in wb:
        for row in ws.iter_rows():
            for cell in row:
                cell.font = font
    wb.save(mergefile)
    
def main():
    print('''
    Script for Merging Data after Docking.

1. Merge the basic information of all crystals in automatedMD/lib/info
2. Merge all docking results data in automatedMD/lib/result
3. Merge all docking with exogenous ligand results data in automatedMD/lib/result

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
    else:
        sys.exit()
    
if __name__ == '__main__':
    main()