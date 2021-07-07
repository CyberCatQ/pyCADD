#!/usr/bin/python
'''

MD轨迹任意指定原子键长键角二面角分析与导出 
Version 1.01

Author YH. W
Last update: 2021/05/21

'''

import pytraj as pt
import xlsxwriter
import os
import re
import sys

try:
    os.makedirs('./pynalysis/dis_ang')
except FileExistsError:
    pass

# Amber Mask
residues = re.compile(r':\d+@\w+')
atoms = re.compile(r'@\w+')
res_num = re.compile(r'(?<=:)\d+(?=@)')
atom_num = re.compile(r'(?<=@)\d+')


def check_mask(lis):  
    '''
    检查Amber Mask是否不合法 

    Paramters
    ----------
    lis: 要检查的AMBER Mask列表

    Return
    ---------
    全部合法返回False 
    否则返回不合法Mask字符串
    '''
    for mask in lis:
        if not re.match(residues, mask) and not re.match(atoms, mask):
            return mask
    return False


def get_mask(flag):
    '''
    获取用户指定的原子AMBER Mask并检查合法性

    Paramters
    ----------
    flag: 判断需要输入Amber Mask数量的标志

    Return
    ----------
    包含全部已输入Amber Mask的列表
    '''

    print('\nAmber Mask Example: :123@CA or @123\n')
    while True:
        mask1 = input('请输入第一个原子Amber Mask:').upper().strip()
        mask2 = input('请输入第二个原子Amber Mask:').upper().strip()
        lis = [mask1, mask2]
        if flag >= 2:
            mask3 = input('请输入第三个原子Amber Mask:').upper().strip()
            lis.append(mask3)
            if flag == 3:
                mask4 = input('请输入第四个原子Amber Mask:').upper().strip()
                lis.append(mask4)
        notpass = check_mask(lis)
        if not notpass:
            break
        else:
            print(notpass, '为不合法的Amber Mask，请重新输入\n')
    return lis


def judge(mask, topfile):
    '''
    判断Amber Mask所属类型： Residue | Atom 表达式
    并从拓扑文件中查询该Amber对象信息 打印

    Paramters
    ----------
    mask: Amber Mask字符串
    topfile: 搜寻用拓扑文件

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


def main(traj, top):
    '''
    用户选择计算模式 调用pytraj计算
    调用xlsxwriter保存数据文件

    Paramters
    ----------
    traj: MD轨迹 Pytraj对象
    top: MD拓扑文件PATH
    '''

    print('''

    请输入计算模式:
    
    1.计算原子间距离(键长)变化
    2.计算键角变化
    3.计算二面角变化

    0.退出

    ''')
    while True:
        flag = input().strip()
        if re.match('[0123]', flag):
            flag = int(flag)
            break
        else:
            print('输入代号无效，请重新输入\n')

    if flag == 0:
        sys.exit()

    lis = get_mask(flag)
    topfile = pt.load_topology(top)
    masks = ''
    for i in lis:
        masks += i + ' '
        judge(i, topfile)

    print('\nProcessing and Saving Data...')
    disangxlsx = xlsxwriter.Workbook(
        './pynalysis/dis_ang/distance_angle_data.xlsx')

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
