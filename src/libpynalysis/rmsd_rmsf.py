#!/usr/bin/python
'''

计算MD轨迹RMSD/F 
Version 1.02

Author YH. W
Last update: 2021/05/21

'''

import pytraj as pt
import os
root_path = os.path.abspath(os.path.dirname(__file__)).split('src')[0]  # 项目路径 绝对路径
pdb_path = root_path.split('automatedMD')[0]                            # PDB项目绝对路径(如果有)
rmsd_path = pdb_path + '/pynalysis/rmsd_rmsf/'

if not os.path.exists(rmsd_path):
    os.makedirs(rmsd_path)

def main(traj, pdb):
    '''
    计算RMSD/RMSF

    Parameters
    ----------
    traj : <object pytraj.traj>
        MD轨迹 Pytraj对象
    pdb : str
        PDB ID字符串
    '''
    
    print('\nStart RMSD/RMSF Analysis...')
    data_rmsd_first = pt.rmsd(traj, ref=0, mask="@CA")  # 计算αC原子RMSD Ref=0以第一帧为参考

    command = "awk '$4 ~ /MOL/ {print $5}' " + pdb_path + pdb + 'comsolvate.pdb'  # 提取位于残基后的配体分子序号
    temp = os.popen(command).readlines()[-1]
    mol_resnum = str(temp).strip()
    res_mask = ':1-' + mol_resnum + '@CA'

    data_rmsf = pt.rmsf(traj, mask=res_mask, options='byres')  # 计算RMSF 所有非水分子

    # 保存数据
    print('RMSD Data:')
    print(data_rmsd_first)
    with open(rmsd_path + 'rmsd.dat', 'w') as rmsdfile:
        for rmsd in data_rmsd_first:
            rmsdfile.write(str(rmsd) + '\n')
        rmsdfile.close()
    print('rmsd.dat saved.', end='\n')
    print('RMSF Data:')
    print(data_rmsf)

    with open(rmsd_path + 'rmsf.dat', 'w') as rmsffile:
        for rmsfdata in data_rmsf:
            rmsffile.write(str(rmsfdata[1]) + '\n')
        rmsffile.close()
    print('rmsf.dat saved.', end='\n')
    print('\n\nRMSD/RMSF Analysis Complete.')
    print(''.center(80,'*'))


if __name__ == '__main__':

    print('请从主程序pynalysis.py中运行此脚本')
