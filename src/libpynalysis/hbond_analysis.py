#!/usr/bin/python
'''

MD轨迹氢键数据分析与导出 
Version 1.10

Author YH. W
Last update: 2021/07/12

'''

import pytraj as pt
import xlsxwriter
import os
root_path = os.path.abspath(os.path.dirname(__file__)).split('src')[0]  # 项目路径 绝对路径
pdb_path = root_path.split('automatedMD')[0]                            # PDB项目绝对路径(如果有)

try:
    os.makedirs(pdb_path + 'pynalysis/hbond')
except FileExistsError:
    pass

def main(traj):
    '''
    氢键分析与数据导出

    Parameters
    ----------
    traj : <object pytraj.traj>
        pytraj MD轨迹对象
    '''
    
    print('\nStart H-Bond Analysis...')
    # hbond 分析
    hb = pt.hbond(traj, mask=':*', distance=4,
                  options='avgout %spynalysis/hbond/avg-hbd.dat printatomnum nointramol' % pdb_path)
    '''
    options可填参数与cpptraj一致
    输出氢键平均信息文件 
    打印原子序号 
    仅计算分子间氢键 
    distance: 识别氢键cutoff值 默认值3埃
    
    '''

    distance_mask = hb.get_amber_mask()[0]
    print('Hbond Distance Mask: {} \n '.format(distance_mask))
    angle_mask = hb.get_amber_mask()[1]
    print('Hbond Hngle Mask: {} \n'.format(angle_mask))

    print("\nHbond Data")
    print(hb.data)  # 1: have hbond; 0: does not have hbond

    # 创建excel工作表
    hbondxlsx = xlsxwriter.Workbook(pdb_path + 'pynalysis/hbond/hbond.xlsx')
    distance_sheet = hbondxlsx.add_worksheet('Distance')
    angle_sheet = hbondxlsx.add_worksheet('Angle')

    # 写入header
    col = 0
    row = 0
    for hbond_mask in distance_mask:
        distance_sheet.write(row, col, hbond_mask)
        col += 1

    col = 0
    row = 0
    for hbond_mask in angle_mask:
        angle_sheet.write(row, col, hbond_mask)
        col += 1

    # 计算氢键距离
    dist = pt.distance(traj, distance_mask)
    print('\nAll Hbond Distance: \n', dist, end='\n')
    # 计算氢键角度
    angle = pt.angle(traj, angle_mask)
    print('\nAll Hbond Angle: \n', angle, end='\n')

    print('\nSaving Data...', end='\n')
    # 写入氢键长度角度数据
    row = 1
    col = 0
    for i in dist[:]:
        for j in i[:]:
            distance_sheet.write(row, col, j)
            row += 1
        col += 1
        row = 1

    row = 1
    col = 0
    for l in angle[:]:
        for k in l[:]:
            angle_sheet.write(row, col, k)
            row += 1
        col += 1
        row = 1

    hbondxlsx.close()

    print('\nDistance and Angle Raw Data(hbond.xlsx) saved.', end='\n')
    print('\nH Bond Analysis Complete\n')
    print(''.center(80,'*'))


if __name__ == '__main__':

    print('请从主程序pynalysis.py中运行此脚本')
