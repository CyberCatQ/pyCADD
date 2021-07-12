#!/usr/bin/python
'''

MD轨迹任意指定原子键长键角二面角分析与导出 
Version 1.10

Author YH. W
Last update: 2021/07/12

'''

import pytraj as pt
import xlsxwriter
import os
import re
import sys
root_path = os.path.abspath(os.path.dirname(__file__)).split('src')[0]  # 项目路径 绝对路径
pdb_path = root_path.split('automatedMD')[0]                            # PDB项目绝对路径(如果有)

try:
    os.makedirs(pdb_path + 'pynalysis/dis_ang')
except FileExistsError:
    pass

# AMBER Mask
residues = re.compile(r':\d+@\w+')
atoms = re.compile(r'@\w+')
res_num = re.compile(r'(?<=:)\d+(?=@)')
atom_num = re.compile(r'(?<=@)\d+')


def check_mask(lis):  
    '''
    检查AMBER Mask是否[不]合法 

    Parameters
    ----------
    lis : list
        要检查的AMBER Mask列表

    Return
    ---------
    False | Invalid String[str]
        全部合法返回False(无不合法字符串)
        否则返回不合法Mask字符串
    '''
    for mask in lis:
        if not re.match(residues, mask) and not re.match(atoms, mask):
            return mask
    return False


def get_mask(flag):
    '''
    获取用户指定的原子AMBER Mask并检查合法性

    Parameters
    ----------
    flag : str
        判断需要输入AMBER Mask数量的标志

    Return
    ----------
    list
        包含全部已输入AMBER Mask的列表
    '''

    print('\nAMBER Mask Example: :123@CA or @123\n')
    while True:
        mask1 = input('Enter the first ATOM AMBER Mask:').upper().strip()
        mask2 = input('Enter the second ATOM AMBER Mask:').upper().strip()
        lis = [mask1, mask2]
        if flag >= 2:
            mask3 = input('Enter the third ATOM AMBER Mask:').upper().strip()
            lis.append(mask3)
            if flag == 3:
                mask4 = input('Enter the forth ATOM AMBER Mask:').upper().strip()
                lis.append(mask4)
        notpass = check_mask(lis)
        if not notpass:
            break
        else:
            print('Invalid AMBER Mask: %s\nPlease Try Again.' % notpass)
    return lis


def judge(mask, topfile):
    '''
    判断AMBER Mask所属类型 ： Residue | Atom 表达式
    并从拓扑文件中查询该AMBER对象信息

    Parameters
    ----------
    mask : str
        AMBER Mask字符串
    topfile: <object pytraj.topology>
        搜寻用拓扑文件对象

    '''
    if re.search(res_num, mask):  # 如果Mask是一个Residue表达式
        resid = int(re.search(res_num, mask).group())-1
        res = topfile.residue(resid)
        atom_index = int(topfile.select(mask)[0])
    elif re.search(atom_num, mask):  # 如果Mask是一个Atom表达式
        atom_index = int(re.search(atom_num, mask).group())-1
        atom = topfile.atom(atom_index)
        res = topfile.residue(int(atom.resid))  # 原子所在分子号与resid相差1
    print('\nAtom ' + mask + ' Info:')
    print('Source Residue:', res.name, res.original_resid)
    print('Atom Index:', atom_index, '\nAtom Type:',
          topfile.atom(atom_index).name, end='\n')


def main(traj, topfile):
    '''
    选择计算模式 调用pytraj计算
    调用xlsxwriter保存数据文件

    Parameters
    ----------
    traj : <object pytraj.traj>
        Pytraj轨迹对象
    topfile : str
        MD拓扑文件PATH
    '''

    print('''

Enter the Calculate Mode Code:
    
1. Calculate the Change in Distance (Bond Length) Between Atoms
2. Calculate the Bond Angle Change
3. Calculate the Dihedral Angle Change

0. Exit

    ''')
    while True:
        flag = input().strip()
        if re.match('[0123]', flag):
            flag = int(flag)
            break
        else:
            print('Invalid Code. Please Try Again.\n')

    if flag == 0:
        sys.exit()

    lis = get_mask(flag)
    topfile = pt.load_topology(topfile)
    masks = ''
    for i in lis:
        masks += i + ' '
        judge(i, topfile)

    print('\nProcessing and Saving Data...')
    disangxlsx = xlsxwriter.Workbook(
        pdb_path + 'pynalysis/dis_ang/distance_angle_data.xlsx')

    def save(sheet, mask, data):  # 保存数据
        sheet.write(0, 0, mask)
        col = 0
        row = 1
        for i in data:
            sheet.write(row, col, i)
            row += 1

    if flag == 1:
        dis = pt.distance(traj, masks)
        print('\nDistance Data:\n', dis)
        distance_sheet = disangxlsx.add_worksheet('Distance')
        save(distance_sheet, masks, dis)

    if flag == 2:
        angle = pt.angle(traj, masks)
        print('Angel Data:\n', angle)
        angel_sheet = disangxlsx.add_worksheet('Angel')
        save(angel_sheet, masks, angle)

    if flag == 3:
        dihedral = pt.dihedral(traj, masks)
        print('Dihedral Data:\n', dihedral)
        dihedral_sheet = disangxlsx.add_worksheet('Dihedral')
        save(dihedral_sheet, masks, dihedral)

    disangxlsx.close()
    print('\ndistance_angle_data.xlsx Data Saved.', end='\n')


if __name__ == '__main__':

    print('请从主程序pynalysis.py中启动此脚本')
