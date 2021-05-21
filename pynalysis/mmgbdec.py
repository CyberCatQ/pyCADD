#!/usr/bin/python
'''

调用MMPBSA.py计算吉布斯自由能变(MM-GB/SA方法) 
Version 1.02

Author YH. W
Last update: 2021/05/21

'''

import os


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
mkdir('./pynalysis/MMPBSA/2gbdec/')


def mmgbdec_prepare():
    '''
    编写MM-GB/SA输入文件
    '''

    text = ''
    text += '#input file for GB and PB calculation\n'
    text += '&general\n'
    text += 'startframe=40010, endframe=50000, interval=10,\n'  # 从40010帧开始至50000帧 间隔为10
    text += 'verbose=2, keep_files=1,\n'  # 打印所有配体 受体 复合物及键合项数据 保留重要临时文件
    text += '/\n'
    text += '&gb\n'
    text += 'igb=5, saltcon=0.15,\n'  # Generalized Born 方法5 | 盐的物质的量浓度0.15M
    text += '/\n'
    text += '&decomp\n'
    text += 'idecomp=1,\n'  # 残基能量分解参数
    text += '/\n'
    with open('./pynalysis/MMPBSA/2gbdec/mmgbdec.in', 'w') as f:
        f.write(text)


def main(mpi_num, pdbcom_prmtop, pdbpro_prmtop, pdblig_prmtop, pdbcomsolvate_prmtop):
    '''
    调用MMPBSA.py.MPI启动MM-GB/SA分析

    Parameters
    ----------
    mpi_num: 使用的核心数量
    pdbcom_prmtop: 复合物拓扑文件PATH
    pdbpro_prmtop: 蛋白拓扑文件PATH
    pdblig_prmtop: 配体拓扑文件PATH
    pdbcomsolvate_prmtop: 溶剂化复合物拓扑文件PATH

    '''
    mmgbdec_prepare()
    print('Start MMGBSA Analysis...')

    command = 'nohup mpirun -np ' + mpi_num + \
        ' MMPBSA.py.MPI -O -do ./pynalysis/MMPBSA/2gbdec/FINAL_DECOMP_MMPBSA.dat ' + \
        '-o ./pynalysis/MMPBSA/2gbdec/FINAL_RESULTS_MMPBSA.dat -i ./pynalysis/MMPBSA/2gbdec/mmgbdec.in' + \
        ' -cp ' + pdbcom_prmtop + ' -rp ' + pdbpro_prmtop + ' -lp ' + pdblig_prmtop + ' -sp ' + \
        pdbcomsolvate_prmtop + ' -y ./npt/npt.mdcrd > ./pynalysis/MMPBSA/2gbdec/mmpbsa.log 2>&1 &'
    os.system(command)

    print('\nMMPBSA.py.MPI jobs submitted.')


if __name__ == '__main__':

    print('请从主程序pynalysis.py中运行此脚本')
