#!/usr/bin/python

import os
import sys
import re
import getopt

root_path = os.path.abspath(os.path.dirname(__file__)).split('src')[0]  # 项目路径 绝对路径
lib_path = root_path + 'lib' + os.sep                                   # 库文件夹路径
doc_path = root_path + 'doc' + os.sep                                   # 文档文件夹路径
pdb_path = root_path.split('automatedMD')[0]                            # PDB项目绝对路径(如果有)
pdb_name = pdb_path.split(os.sep)[-2]                                   # PDB项目名称(如果有)
src_path = os.path.dirname(__file__)                                    # 源代码文件夹路径

class MDprepare:
    '''

Python MD Prepare 主程序 
Version 1.00

Author YH. W
Last Update: 2021/07/09

'''
    def __init__(self, pdbid) -> None:
        self.pdbid = pdbid          # PDB ID
        self.pro4leap_file = ''     # 准备好的LEaP可读蛋白文件
        self.lig4leap_file = ''     # 准备好的LEaP可读配体文件
        self.mpi_num = 4            # 使用的CPU核心数
        self.ele = 0                # 配体分子电荷数
        self.spin_multiplicity = 1  # 配体分子自旋多重度

    def get_pdbid(self):
        '''
        获取用户输入的PDBID并检查合法性

        Return
        ----------
        str
            PDB ID字符串

        '''
        def check(pdb):
            return re.fullmatch(r'^\d[0-9a-zA-Z]{3,}$', pdb)

        if self.pdbid:
            pdb = self.pdbid
        else:
            pdb = pdb_name

        if check(pdb):
            self.pdbid = pdb
            return pdb
        else:
            while True:
                pdb = str(
                    input('\nTo get PDB ID automatically, please change the name of crystal folder to PDBID\n Input PDB ID:')).strip().upper()
                if check(pdb):
                    self.pdbid = pdb
                    return pdb
                else:
                    print('Please enter the correct PDBID. Try Again.')

    def protein_prepare(self, pdbid=None):
        '''
        受体蛋白文件预处理 : \n
        PDB文件格式化 for Amber | 去除原生H原子 | 使用rudece添加H原子 | 再次格式化

        Temporary files:
            PDBpro_dry.pdb 
                原始蛋白去除水分子 添加缺失重原子\n      
            PDBpro_noH.pdb
                pro_dry.pdb结构继续去除所有H原子\n     
            PDBpro_leap.pdb
                pro_noH.pdb结构格式化为LEaP可识别格式\n      

        Parameters
        ----------
        pdb : str
            PDB ID字符串

        '''
        if not pdbid:
            pdbid = self.pdbid
        if not os.path.exists('%spro.pdb' % pdbid):
            raise FileNotFoundError('%spro.pdb is not found in %s' % (pdbid, pdb_path))
        os.system(
            'pdb4amber -i %spro.pdb -o %spro_dry.pdb -p --dry --add-missing-atoms ' % (pdbid, pdbid))
        os.system('pdb4amber -i %spro_dry.pdb -o %spro_noH.pdb -y' %
                (pdbid, pdbid))  # 必须分两步否则H原子删除不完全
        os.system('pdb4amber -i %spro_noH.pdb -o %spro_leap.pdb' % (pdbid, pdbid))

        pro4leap_file = '%spro_leap.pdb' % pdbid
        if os.path.exists(pro4leap_file):
            self.pro4leap_file = pro4leap_file
        else:
            raise RuntimeError('Protein Processing Failed.')

        print('\nProtein Process Complete. %s savad.' % pro4leap_file)
        print(''.center(80, '*'), end='\n')


    def ligand_prepare(self, pdbid=None, mpi_num=4, ele=0, spin_multiplicity=1):
        '''
        小分子文件准备:高斯坐标优化与RESP2(0.5)电荷计算

        Parameters
        ----------
        pdb : str 
            PDB ID字符串或配体小分子名
        mpi_num : int | str
            需要使用的CPU核心数 默认为4
        ele : int | str
            分子电荷数 默认为0
        spin_multiplicity : int | str
            自旋多重度 默认为1

        '''
        if not pdbid:
            pdbid = self.pdbid

        print('\nPrepare to Runnning Gaussian...\n')
        # 准备计算RESP电荷 调用Gaussian与Multiwnf
        gauss_path = os.environ['g16root']

        with open(gauss_path+'/g16/Default.Route', 'w') as gaussain:  # 修改高斯输入文件头 w模式覆盖原文件
            gaussain.write('-M- 16GB\n')                              # 调用的内存大小 默认16GB
            gaussain.write('-P- %s ' % mpi_num)                       # 使用的CPU核心数 默认4

        print('\nPrepare to Calculate...\n')
        os.system('chmod 777 ./RESP2.sh && ./RESP2.sh %slig.pdb %s %s' % (pdbid, ele,
                spin_multiplicity))  # 输入参数调用Multiwnf计算RESP2电荷 自动输出为当前目录的同名lig.pqr文件 包含高斯优化坐标&RRESP2电荷

        # pqr文件无法继续进行prepin文件生成 需要通过openbabel转换为mol2
        print('\nPrepare to Generate Parmter File...\n')

        os.system('obabel -ipqr %slig.pqr -omol2 -O %slig.mol2' % (pdbid, pdbid))   # 依赖openbabel软件
        if not os.path.exists('%slig.mol2' % pdbid):
            raise RuntimeError('Openbabel is not installed. Please install it using the following command:\nsudo apt install openbabel')

        os.system('antechamber -fi mol2 -i %slig.mol2 -fo prepi -o %slig.prepin' %
                (pdbid, pdbid))  # antechamber转换为Amber prepi文件
        os.system(
            'antechamber -fi mol2 -i %slig.mol2 -o %slig_leap.pdb -fo pdb' % (pdbid, pdbid))    # 生成LEaP可读配体小分子的PDB文件
        os.system('parmchk2 -i %slig.prepin -f prepi -o %slig.frcmod' %
                (pdbid, pdbid))  # 生成Amber Parmter 参数文件

        prepin_file = '%slig.prepin' % pdbid
        frcmod_file = '%slig.frcmod' % pdbid
        lig4leap_file = '%slig_leap.pdb' % pdbid
        if not os.path.exists(prepin_file) or not os.path.exists(frcmod_file) or not os.path.exists(lig4leap_file):
            raise RuntimeError('Antechamber module running failed or AmberTools is not installed.')
        
        self.lig4leap_file = lig4leap_file
        print('\nLigand Preparation Progress Complete.')
        print(''.center(80, '*'), end='\n')


    def leap_prepare(self, pdbid=None):
        '''
        创建LEaP输入文件

        Parameters
        ----------
        pdb : str
            PDB ID或小分子名

        '''
        if not pdbid:
            pdb = self.pdbid
        else:
            pdb = pdbid

        # MD核心参数 请参照AMBER SANDER文档进行调试
        with open('tleaplig.in', 'w') as f:

            text = 'source leaprc.protein.ff14SB\n'     # 蛋白力场
            text += 'source leaprc.gaff2\n'             # 
            text += 'source leaprc.water.tip3p\n'       # 水箱模型

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

    def __usage():
        '''
Preprocessing for Molecular Dynamics Simulation.

Usage: 
pyMDprepare.py -n <number of cpu> [ -e <charge> -s <spin multiplicity> ]

Descrption:
        
    -h      HELP        
    -n      Number of CPU cores to use      
    -e      Ligand molecular charge number      (Defualt: 0)
    -s      Ligand molecular spin multiplicity  (Default: 1)  
'''

    def __procsee_argvs(self):
        '''
        解析命令行参数
        '''
        try:
            opts, argvs = getopt.getopt(sys.argv[1:], '-he:s:n:')
        except getopt.GetoptError as err:
            print(str(err))
            print(self.__usage.__doc__)
            sys.exit(1)

        for opt, arg in opts:

            if opt == '-h':
                print(self.__usage.__doc__)
                sys.exit(1)
            elif opt == '-n':
                self.mpi_num = arg
            elif opt == '-e':
                self.ele = arg
            elif opt == '-s':
                self.spin_multiplicity = arg

    def __check_requirements():
        '''
        检查必要文件与依赖
        '''
        if not os.path.exists(src_path + '/RESP2.sh'):  # 检查必要文件与依赖程序
            raise RuntimeError('Script RESP2.sh is not Found!')

        try:
            Multiwfn = os.environ['Multiwfnpath']       # 检查Multiwfn
        except KeyError:
            raise RuntimeError('Multiwfn is not installed.')

        try:                                            # 检查高斯软件路径
            gauss_path = os.environ['g16root']
        except KeyError:                      
            while True:
                print('Can not find directory of Gaussin Software.\nPlease specify g16root directory: ')
                gauss_path = os.path.abspath(input())
                if not os.path.exists(gauss_path + '/g16'):
                    print('Please enter the correct directory path. Try again.')
                else:
                    break

    def main(self):
        '''
        执行MD Prepare流程
        '''
        pdb = self.get_pdbid()
        self.__check_requirements()
        self.__procsee_argvs()

        mpi_num = self.mpi_num
        ele = self.ele
        spin_multiplicity = self.spin_multiplicity

        print('\n\n')
        print('Python Script For Amber MD Prepatation'.center(
            80, '-'), end='\n')  # Title
        print('\n Processing Entry ID: %s' % pdb, end='\n')

        self.protein_prepare(pdb)
        self.ligand_prepare(pdb, mpi_num, ele, spin_multiplicity)
        self.leap_prepare(pdb)

        print('\nRunning tleap...\n')
        os.system('tleap -f tleaplig.in')

        os.system('ambpdb -p %scomsolvate.prmtop < %scomsolvate.inpcrd > %scomsolvate.pdb' % (
            pdb, pdb, pdb))  # 生成溶剂化的复合物PDB文件

        print('\nMD Preparation Process Complete.\n')
        print(''.center(80, '-'), end='\n')


if __name__ == '__main__':
    console = MDprepare()
    if len(sys.argv) == 1:
        print(console.__usage.__doc__)
    else:
        console.main()
