#!/usr/bin/python
'''

Python分析MD轨迹 主程序 version 1.02

Author YH. W
Last Update: 2021/05/21

'''

import os
import re
import sys
import multiprocessing
cwd = str(os.getcwd())
sys.path.append(cwd)

try:
    from pynalysis import rmsd_rmsf, hbond_analysis, lowest_energy_structure, mmgbdec, nmode, dis_ang

except ImportError:
    print('Error: 请将py脚本文件与pynalysis包置于同一文件夹内')
    sys.exit(2)

try:
    import xlsxwriter
    import pytraj as pt
except ImportError:
    print('''\n
    Error: 缺失依赖包
    请使用pip install安装依赖包后重试\n
    ''')
    raise


def get_pdb():
    '''
    尝试自动获取PDB ID 
    否则提示用户手动输入并检查合法性

    Return
    ----------
    PDB ID字符串
    '''
    def check(pdb):
        return re.fullmatch(r'^\d[0-9a-zA-Z]{3,}$', pdb)

    pdb = str(os.path.split(cwd)[1])
    if check(pdb):
        return pdb
    else:
        while True:
            pdb = str(
                input('\n要自动获取PDB代号 请将晶体文件夹修改为晶体ID\n请手动输入晶体PDB代号:')).strip().upper()
            if check(pdb):
                return pdb
            else:
                print('请输入正确的PDB代号!')


def load_traj(top):
    '''
    从拓扑文件读取信息
    并调用pytraj加载轨迹文件
    轨迹文件默认PATH ./npt/npt.mdcrd

    Parameters
    ----------
    top: MD轨迹拓扑文件PATH

    Return
    ----------
    pytraj MD轨迹对象
    '''

    print('\nStart Loading Trajectory...', end='\n')
    print('\nLoading Trajectory File: ./npt/npt.mdcrd',
          end='\n')  # 轨迹文件默认为npt.mdcrd
    traj = pt.iterload('./npt/npt.mdcrd', top)  # 读取轨迹文件与拓扑文件
    print('\nTrajectory Loading Complete.\n \nBasic info:', end='\n')
    print(traj, end='\n')  # 显示基本信息

    return traj


def run():
    '''
    检查必需文件并启动MD轨迹分析
    '''

    pdb = get_pdb()
    pdbcom_prmtop = pdb + 'com.prmtop'
    pdbpro_prmtop = pdb + 'pro.prmtop'
    pdblig_prmtop = pdb + 'lig.prmtop'
    pdbcomsolvate_prmtop = pdb + 'comsolvate.prmtop'
    traj_file = './npt/npt.mdcrd'

    # 检查所需文件是否存在
    try:
        print('\nChecking Needed Files...')
        for i in [pdbcom_prmtop, pdbpro_prmtop, pdblig_prmtop, pdbcomsolvate_prmtop, traj_file]:
            f = open(i)
            f.close()
    except IOError:
        print('Error: 无法读取文件 %s 或文件不存在' % i)
        sys.exit(2)

    else:
        dic = dict(PDB代号=pdb, 复合物拓扑文件=pdbcom_prmtop, 蛋白拓扑文件=pdblig_prmtop,
                   配体拓扑文件=pdblig_prmtop, 溶剂化复合物拓扑文件=pdbcomsolvate_prmtop, 轨迹文件=traj_file)
        key = ['PDB代号', '复合物拓扑文件', '蛋白拓扑文件', '配体拓扑文件', '轨迹文件']
        for k in key:
            print('%s:%s' % (k, dic[k]), end='\n')

        print('\nAll Needed Files Dectected.')

    traj = load_traj(pdbcomsolvate_prmtop)
    flag = input('''
    请输入需要进行的分析操作代号:

    1.RMSD & RMSF计算
    2.H-Bond 距离与角度分析
    3.最低能量构象提取
    4.MMPBSA.py并行计算 吉布斯自由能变
    5.MMPBSA.py并行计算 熵变
    6.指定原子键长、键角、二面角变化分析
    
    0.退出

    ''')

    record = []

    if '0' in flag:
        sys.exit(0)
    if '1' in flag:
        process1 = multiprocessing.Process(
            target=rmsd_rmsf.main, args=(traj, pdb))
        process1.start()
        record.append(process1)
    if '2' in flag:
        process2 = multiprocessing.Process(
            target=hbond_analysis.main, args=(traj,))
        process2.start()
        record.append(process2)
    if '3' in flag:
        process3 = multiprocessing.Process(
            target=lowest_energy_structure.main, args=(traj,))
        process3.start()
        record.append(process3)

    if '4' in flag or '5' in flag:
        mpi_num = input('输入MMGBSA进程要使用的CPU核心数:')
    if '4' in flag:
        process4 = multiprocessing.Process(target=mmgbdec.main, args=(mpi_num,pdbcom_prmtop, pdbpro_prmtop, pdblig_prmtop, pdbcomsolvate_prmtop))
        process4.start()
        record.append(process4)
    if '5' in flag:
        process5 = multiprocessing.Process(target=nmode.main, args=(mpi_num,pdbcom_prmtop, pdbpro_prmtop, pdblig_prmtop, pdbcomsolvate_prmtop))
        process5.start()
        record.append(process5)

    elif flag == '6':   #需要用户手动操作 无法多进程化
        dis_ang.main(traj, pdbcomsolvate_prmtop)
    else:
        print('Error: 输入代号无效')
        sys.exit(2)

    for process in record:
        process.join()

    print('\nAll Analysis Processing Complete.')
    print(''.center(80, '-'))


if __name__ == '__main__':
    run()
