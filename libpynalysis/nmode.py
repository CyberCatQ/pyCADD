#!/usr/bin/python
'''
调用MMPBSA.py计算熵变(nmode方法) 
Version 1.10

Author YH. W
Last update: 2021/07/12

'''

import os
root_path = os.path.abspath(os.path.dirname(__file__)).split('src')[0]  # 项目路径 绝对路径
pdb_path = root_path.split('automatedMD')[0]                            # PDB项目绝对路径(如果有)
MMPBSA_path = pdb_path + 'pynalysis/MMPBSA/'

def mkdir(path):
    '''
    创建必要文件夹

    Parameters
    ----------
    path : str
        文件夹PATH
    '''

    try:
        os.makedirs(path)
    except FileExistsError:
        pass


mkdir(MMPBSA_path)
mkdir(MMPBSA_path + 'nmode/')


def nmode_prepare():
    '''
    编写MM-GB/SA熵计算输入文件
    '''

    text = ''
    text += '#input file with nmode analysis\n'
    text += '&general\n'
    text += 'startframe=40100, endframe=50000, interval=100,\n'  # 从40010帧开始至50000帧 间隔为100
    text += 'verbose=2, keep_files=1,\n'  # 打印所有配体 受体 复合物及键合项数据 保留重要临时文件
    text += '/\n'
    text += '&nmode\n'
    # 能量最小化循环最多10000 | 能量最小化收敛准则0.001 | Generalized Born model:HCT GB | 离子强度0.15M
    text += 'maxcyc=10000, drms=0.001, nmode_igb=1, nmode_istrng=0.15,\n'
    text += '/\n'

    with open(MMPBSA_path + 'nmode/nmode.in', 'w') as f:
        f.write(text)


def main(mpi_num, pdbcom_prmtop, pdbpro_prmtop, pdblig_prmtop, pdbcomsolvate_prmtop):
    '''
    调用MMPBSA.py.MPI启动熵计算分析

    Parameters
    ----------
    mpi_num : int | str
        使用的核心数量
    pdbcom_prmtop : str
        复合物拓扑文件PATH
    pdbpro_prmtop : str
        蛋白拓扑文件PATH
    pdblig_prmtop : str
        配体拓扑文件PATH
    pdbcomsolvate_prmtop : str
        溶剂化复合物拓扑文件PATH

    '''

    nmode_prepare()
    print('Start MMGBSA(Entropy) Analysis...')

    command = 'nohup mpirun -np ' + mpi_num + ' MMPBSA.py.MPI -O -o ' + MMPBSA_path + 'nmode/FINAL_RESULTS_MMPBSA.dat -i ' + MMPBSA_path + 'nmode/nmode.in' + ' -cp ' + pdbcom_prmtop + \
        ' -rp ' + pdbpro_prmtop + ' -lp ' + pdblig_prmtop + ' -sp ' + pdbcomsolvate_prmtop + \
        ' -y ' + pdb_path + 'npt/npt.mdcrd > ' + MMPBSA_path + 'nmode/nmode.log 2>&1 &'
    os.system(command)

    print('MMPBSA.py.MPI jobs submitted.')


if __name__ == '__main__':

    print('请从主程序pynalysis.py中运行此脚本')
