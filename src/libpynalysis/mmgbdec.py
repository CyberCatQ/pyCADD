#!/usr/bin/python
'''

调用MMPBSA.py计算吉布斯自由能变(MM-GB/SA方法) 
Version 1.10

Author YH. W
Last update: 2021/07/12

'''

import os

root_path = os.path.abspath(os.path.dirname(__file__)).split('src')[0]  # 项目路径 绝对路径
pdb_path = root_path.split('automatedMD')[0]                            # PDB项目绝对路径(如果有)

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


mkdir('%spynalysis/MMPBSA/' % pdb_path)
mkdir('%spynalysis/MMPBSA/gbdec/' % pdb_path)


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
    with open('%spynalysis/MMPBSA/gbdec/mmgbdec.in' % pdb_path, 'w') as f:
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

    command = 'nohup mpirun -np %s MMPBSA.py.MPI -O -do %spynalysis/MMPBSA/gbdec/FINAL_DECOMP_MMPBSA.dat -o %spynalysis/MMPBSA/gbdec/FINAL_RESULTS_MMPBSA.dat -i %spynalysis/MMPBSA/gbdec/mmgbdec.in -cp %s -rp %s -lp %s -sp %s -y %snpt/npt.mdcrd > %spynalysis/MMPBSA/gbdec/mmpbsa.log 2>&1 &' % (
        mpi_num, pdb_path, pdb_path, pdb_path, pdbcom_prmtop, pdbpro_prmtop, pdblig_prmtop, pdbcomsolvate_prmtop, pdb_path, pdb_path)
    os.system(command)

    print('\nMMPBSA.py.MPI jobs submitted.')


if __name__ == '__main__':

    print('请从主程序pynalysis.py中运行此脚本')
