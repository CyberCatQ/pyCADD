#!/usr/bin/python
'''

计算MD轨迹RMSD/F vesion 1.02

Author YH. W
Last update: 2021/05/21

'''

import pytraj as pt
import os
if not os.path.exists('./pynalysis/rmsd_rmsf/'):
    os.makedirs('./pynalysis/rmsd_rmsf/')

def main(traj, pdb):
    '''
    计算RMSD/RMSF

    Parameters
    ----------
    traj: MD轨迹 Pytraj对象
    pdb: PDB ID字符串
    '''
    
    print('\nStart RMSD/RMSF Analysis...')
    data_rmsd_first = pt.rmsd(traj, ref=0, mask="@CA")  # 计算αC原子RMSD Ref=0以第一帧为参考

    command = "awk '$4 ~ /MOL/ {print $5}' ./" + pdb + 'comsolvate.pdb'  # 提取位于残基后的配体分子序号
    temp = os.popen(command)
    mol_resnum = str(temp.readlines()[-1]).strip()
    res_mask = ':1-' + mol_resnum + '@CA'

    data_rmsf = pt.rmsf(traj, mask=res_mask, options='byres')  # 计算RMSF 所有非水分子

    # 保存数据
    print('RMSD Data:')
    print(data_rmsd_first)
    with open('./pynalysis/rmsd_rmsf/rmsd.dat', 'w') as rmsdfile:
        for rmsd in data_rmsd_first[:]:
            rmsdfile.write(str(rmsd) + '\n')
        rmsdfile.close()
    print('rmsd.dat saved.', end='\n')
    print('RMSF Data:')
    print(data_rmsf)

    with open('./pynalysis/rmsd_rmsf/rmsf.dat', 'w') as rmsffile:
        for rmsfdata in data_rmsf[:]:
            rmsffile.write(str(rmsfdata[1]) + '\n')
        rmsffile.close()
    print('rmsf.dat saved.', end='\n')
    print('\n\nRMSD/RMSF Analysis Complete.')
    print(''.center(80,'*'))


if __name__ == '__main__':

    print('请从主程序pynalysis.py中运行此脚本')
