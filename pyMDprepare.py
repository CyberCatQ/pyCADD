#!/usr/bin/python
'''

Python MD Prepare 主程序 
Version 1.00

Author YH. W
Last Update: 2021/05/18

'''

import os
import sys
import re
import getopt

cwd = str(os.getcwd())
sys.path.append(cwd)


def get_pdbid():
    '''
    获取用户输入的PDBID并检查合法性

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
                input('\n要自动获取PDB ID 请将晶体文件夹修改为PDB ID\n请手动输入晶体PDB ID:')).strip().upper()
            if check(pdb):
                return pdb
            else:
                print('请输入正确的PDB ID!')


def protein_prepare(pdb):
    '''
    受体蛋白文件预处理: PDB文件格式化 for Amber | 去除原生H原子 | 使用rudece添加H原子 | 再次格式化

    Parameters
    ----------
    pdb: PDB ID字符串

    '''

    os.system(
        'pdb4amber -i %spro.pdb -o %spro_dry.pdb -p --dry --add-missing-atoms ' % (pdb, pdb))
    os.system('pdb4amber -i %spro_dry.pdb -o %spro_noH.pdb -y' %
              (pdb, pdb))  # 必须分两步否则H原子删除不完全
    os.system('pdb4amber -i %spro_noH.pdb -o %spro_leap.pdb' % (pdb, pdb))

    print('\nProtein Process Complete.')
    print(''.center(80, '*'), end='\n')


def ligand_prepare(pdb, mpi_num, ele, spin_multiplicity):
    '''
    小分子文件准备:高斯坐标优化与RESP2(0.5)电荷计算

    Parameters
    ----------
    pdb: PDB ID字符串或配体小分子名
    mpi_num: 需要使用的CPU核心数
    ele: 分子电荷数
    spin_multiplicity: 自旋多重度

    '''
    print('\nPrepare to Runnning Gaussian...\n')
    # 准备计算RESP电荷 调用Gaussian与Multiwnf
    gauss_path = os.environ['g16root']
    with open(gauss_path+'/g16/Default.Route', 'w') as gaussain:  # 修改高斯输入文件头 w模式覆盖原文件
        gaussain.write('-M- 16GB\n')
        gaussain.write('-P- %s ' % mpi_num)

    print('\nPrepare to Calculate...\n')
    os.system('chmod 777 ./RESP2.sh && ./RESP2.sh %slig.pdb %s %s' % (pdb, ele,
              spin_multiplicity))  # 输入参数调用Multiwnf计算RESP2电荷 自动输出为当前目录的同名lig.pqr文件 包含高斯优化坐标&RRESP2电荷

    # pqr文件无法继续进行prepin文件生成 需要通过openbabel转换为mol2
    print('\nPrepare to Generate Parmter File...\n')

    os.system('obabel -ipqr %slig.pqr -omol2 -O %slig.mol2' % (pdb, pdb))
    os.system('antechamber -fi mol2 -i %slig.mol2 -fo prepi -o %slig.prepin' %
              (pdb, pdb))  # antechamber转换为Amber prepi文件
    os.system(
        'antechamber -fi mol2 -i %slig.mol2 -o %slig_leap.pdb -fo pdb' % (pdb, pdb))
    os.system('parmchk2 -i %slig.prepin -f prepi -o %slig.frcmod' %
              (pdb, pdb))  # 生成Amber Parmter 参数文件

    print('\nLigand Preparation Progress Complete.')
    print(''.center(80, '*'), end='\n')


def leap_prepare(pdb):
    '''
    创建LEaP输入文件

    Parameters
    ----------
    pdb: PDB ID或小分子名

    '''

    with open('tleaplig.in', 'w') as f:

        text = 'source leaprc.protein.ff14SB\n'
        text += 'source leaprc.gaff2\n'
        text += 'source leaprc.water.tip3p\n'

        text += 'Ligand'.center(80, '#')
        text += '\nloadamberprep %slig.prepin\n' % pdb
        text += 'loadamberparams %slig.frcmod\n' % pdb
        text += 'lig = loadpdb %slig_leap.pdb\n' % pdb
        text += 'saveamberparm lig %slig.prmtop %slig.inpcrd\n' % (pdb, pdb)

        text += 'Protein'.center(80, '#')
        text += '\npro = loadpdb %spro_leap.pdb\n' % pdb
        text += 'saveamberparm pro %spro.prmtop %spro.inpcrd\n' % (pdb, pdb)

        text += 'COM = combine { pro lig }\n'
        text += 'savepdb COM %scom.pdb\n' % pdb

        text += 'Complex'.center(80, '#')
        text += '\ncom = loadpdb %scom.pdb\n' % pdb
        text += 'saveamberparm com %scom.prmtop %scom.inpcrd\n' % (pdb, pdb)
        text += 'solvatebox com TIP3PBOX 12.0 0.75\n'   #水箱大小太小将导致MD模拟错误
        text += 'check com\n'
        text += 'addions com Na+ 0\n'
        text += 'addions com Cl- 0\n'
        text += 'saveamberparm com %scomsolvate.prmtop %scomsolvate.inpcrd\n' % (
            pdb, pdb)
        text += ''.center(80, '#')
        text += '\nquit\n'

        f.write(text)

    print('\nLEaP Input File Generated.')
    print(''.center(80, '*'), end='\n')


def main(argv):
    '''
    主函数
    获取命令行参数并执行MD Preparation流程

    Parameters
    ----------
    -h 使用帮助
    -n 需要使用的CPU核心数量
    -e 分子电荷数(Defualt: 0)
    -s 自旋多重度(Default: 1)  

    '''

    mpi_num = ''
    ele = '0'
    spin_multiplicity = '1'

    try:
        opts, argvs = getopt.getopt(argv, '-he:s:n:')
    except getopt.GetoptError:
        print(
            "usage: pyMDprepare.py -n <number of cpu> [ -e <charge> -s <spin multiplicity> ]\nUse pyMDprepare.py -h for more information.")
        sys.exit(2)

    for opt, arg in opts:

        if opt == '-h':
            print(
                'usage: pyMDprepare.py -n <number of cpu> [ -e <charge> -s <spin multiplicity> ]')
            print('''
    Options
    ----------
    -n 需要使用的CPU核心数量
    -e 分子电荷数(Defualt: 0)
    -s 自旋多重度(Default: 1)
            ''')
            sys.exit()

        elif opt == '-n':
            mpi_num = arg
        elif opt == '-e':
            ele = arg
        elif opt == '-s':
            spin_multiplicity = arg

    if not os.path.exists('RESP2.sh'):  # 检查必要文件与依赖程序
        print('Error: Script RESP2.sh is not Found!')
        sys.exit(2)
    try:
        Multiwfn = os.environ['Multiwfnpath']
        g16 = os.environ['g16root']
    except KeyError:
        print('\n缺失依赖程序!\n')
        raise

    pdb = get_pdbid()
    print('\n\n')
    print('Python Script For Amber MD Prepatation'.center(
        80, '-'), end='\n')  # Title
    print('\n Processing Entry ID: %s' % pdb, end='\n')

    protein_prepare(pdb)
    ligand_prepare(pdb, mpi_num, ele, spin_multiplicity)
    leap_prepare(pdb)

    print('\nRunning tleap...\n')
    os.system('tleap -f tleaplig.in')

    os.system('ambpdb -p %scomsolvate.prmtop < %scomsolvate.inpcrd > %scomsolvate.pdb' % (
        pdb, pdb, pdb))  # 生成溶剂化的复合物PDB文件

    print('\nMD Preparation Process Complete.\n')
    print(''.center(80, '-'), end='\n')


if __name__ == '__main__':
    if len(sys.argv) == 1:
        print(
            "\nusage: pyMDprepare.py -n <number of cpu> [ -e <charge> -s <spin multiplicity> ]\n\nUse pyMDprepare.py -h for more information.\n")
        sys.exit()
    else:
        main(sys.argv[1:])
