#!/usr/bin/python
'''
调用MMPBSA.py计算熵变(nmode方法) 
vesion 1.01

Author YH. W
Last update: 2021/05/21

'''

import os

from matplotlib.pyplot import text


def mkdir(path):
    '''
    创建必要文件夹

    Paramters
    ----------
    path: 文件夹PATH
    '''

    try:
        os.makedirs(path)
    except FileExistsError:
        pass


mkdir('./pynalysis/MMPBSA/')
mkdir('./pynalysis/MMPBSA/3nmode/')


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

    with open('./pynalysis/MMPBSA/3nmode/nmode.in', 'w') as f:
        f.write(text)


def main(mpi_num, pdbcom_prmtop, pdbpro_prmtop, pdblig_prmtop, pdbcomsolvate_prmtop):
    '''
    调用MMPBSA.py.MPI启动熵计算分析

    Parameters
    ----------
    mpi_num: 使用的核心数量
    pdbcom_prmtop: 复合物拓扑文件PATH
    pdbpro_prmtop: 蛋白拓扑文件PATH
    pdblig_prmtop: 配体拓扑文件PATH
    pdbcomsolvate_prmtop: 溶剂化复合物拓扑文件PATH

    '''

    nmode_prepare()
    print('Start MMGBSA(Entropy) Analysis...')

    command = 'nohup mpirun -np ' + mpi_num + \
        ' MMPBSA.py.MPI -O -o ./pynalysis/MMPBSA/3nmode/FINAL_RESULTS_MMPBSA.dat -i ./pynalysis/MMPBSA/3nmode/nmode.in' + \
        ' -cp ' + pdbcom_prmtop + ' -rp ' + pdbpro_prmtop + ' -lp ' + pdblig_prmtop + \
        ' -sp ' + pdbcomsolvate_prmtop + ' -y ./npt/npt.mdcrd > ./mmpbsa.log 2>&1 &'
    os.system(command)

    print('MMPBSA.py.MPI jobs submitted.')


if __name__ == '__main__':

    print('请从主程序pynalysis.py中运行此脚本')
