#!/usr/bin/python
import os
import re
import sys
import multiprocessing
import xlsxwriter           # 无依赖包将报错     pip install it
import pytraj as pt         # 无依赖包将报错

cwd = str(os.getcwd())
sys.path.append(cwd)
root_path = os.path.abspath(os.path.dirname(__file__)).split('src')[0]  # 项目路径 绝对路径
lib_path = root_path + 'lib' + os.sep                                   # 库文件夹路径
doc_path = root_path + 'doc' + os.sep                                   # 文档文件夹路径
pdb_path = root_path.split('automatedMD')[0]                            # PDB项目绝对路径(如果有)
pdb_name = os.path.basename(pdb_path)                                   # PDB项目名称(如果有)

try:
    from libpynalysis import rmsd_rmsf, hbond_analysis, lowest_energy_structure, mmgbdec, nmode, dis_ang
except ImportError:
    raise ImportError('Please place this script file in the same folder as the libpynalysis package.')

class Pynalysis:
    '''

Python分析MD轨迹 主程序 
Version 1.10

Author YH. W
Last Update: 2021/07/12

'''

    def __init__(self, pdbid) -> None:
        self.pdbid = pdbid                      # PDBID
        self.top_file = ''
        self.traj_file = './npt/npt.mdcrd'      # 默认轨迹文件PATH
        self.pdbcom_prmtop = ''                 # 复合物拓扑文件
        self.pdbpro_prmtop = ''                 # 蛋白拓扑文件
        self.pdblig_prmtop = ''                 # 配体拓扑文件
        self.pdbcomsolvate_prmtop = ''          # 溶剂化复合物拓扑文件

    @staticmethod
    def __check_file(path):
        '''
        检查文件是否存在

        Parameters
        ----------
        path : str
            文件路径
        '''
        if not os.path.exists(path):
            raise RuntimeError('%s is not found' % path)

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


    def load_traj(self, topfile, trajfile):
        '''
        从拓扑文件读取信息
        并调用pytraj加载轨迹文件

        Parameters
        ----------
        topfile : str
            拓扑文件PATH
        trajfile : str
            轨迹文件PATH

        Return
        ----------
        <object pytraj.traj>
            pytraj MD轨迹对象
        '''
        if not topfile:
            topfile = self.top_file
        if not trajfile:
            trajfile = self.traj_file
        self.__check_file(topfile)
        self.__check_file(trajfile)
        top_abs = os.path.abspath(topfile)
        traj_abs = os.path.abspath(trajfile)

        print('\nStart Loading Trajectory...', end='\n')
        print('Topology File: %s' % top_abs, end='\n')               # 拓扑文件
        print('Trajectory File: %s' % traj_abs, end='\n')            # 轨迹文件
        traj = pt.iterload(traj_abs, top_abs)                        # 读取轨迹文件与拓扑文件
        print('\nTrajectory Loading Complete.\n \nBasic info:', end='\n')
        print(traj, end='\n')  # 显示基本信息

        return traj

    def check_need_file(self):
        '''
        检查所需文件是否存在
        '''
        pdb = self.pdbid
        self.__check_file(self.pdbcom_prmtop)
        self.__check_file(self.pdbpro_prmtop)
        self.__check_file(self.pdblig_prmtop)
        self.__check_file(self.pdbcomsolvate_prmtop)
        self.__check_file(self.traj_file)
        pdbcom_prmtop = self.pdbcom_prmtop
        pdbpro_prmtop = self.pdbpro_prmtop
        pdblig_prmtop = self.pdblig_prmtop
        pdbcomsolvate_prmtop = self.pdbcomsolvate_prmtop
        traj_file = self.traj_file

        dic = dict(PDBID=pdb, Complex_Topology_File=pdbcom_prmtop, Protein_Topology_File=pdbpro_prmtop,
                Ligand_Topology_File=pdblig_prmtop, SolvatedCom_Topology_File=pdbcomsolvate_prmtop, Trajectory_File=traj_file)
        key = ['PDBID', 'Complex_Topology_File', 'Protein_Topology_File', 'Ligand_Topology_File', 'SolvatedCom_Topology_File', 'Trajectory_File']
        for k in key:
            print('%s:%s' % (k, dic[k]), end='\n')

        print('\nAll Needed Files Dectected.')

    def main(self):
        '''
        启动MD轨迹分析
        '''
        os.chdir(pdb_path)
        pdb = self.get_pdbid()

        self.pdbcom_prmtop = pdb + 'com.prmtop'
        self.pdbpro_prmtop = pdb + 'pro.prmtop'
        self.pdblig_prmtop = pdb + 'lig.prmtop'
        self.pdbcomsolvate_prmtop = pdb + 'comsolvate.prmtop'
        self.check_need_file()

        pdbcom_prmtop = self.pdbcom_prmtop
        pdbpro_prmtop = self.pdbpro_prmtop
        pdblig_prmtop = self.pdblig_prmtop
        pdbcomsolvate_prmtop = self.pdbcomsolvate_prmtop
        traj_file = self.traj_file
        traj = self.load_traj(pdbcomsolvate_prmtop, traj_file)
        
        flag = input('''
        Please Enter the Code of Analysis to be performed:

        1. RMSD & RMSF Calculate
        2. H-Bond Distance and Angle Analysis
        3. Extract The Lowest Energy Conformation
        4. MMPBSA.py Calculate in parallel: Enthalpy
        5. MMPBSA.py Calculate in parallel: Entropy
        6. Bond Length, Bond Angle, and Dihedral Angle Change Analysis of Specific Atoms
        
        0. Exit

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
            mpi_num = input('Please Enter the Number of CPU used by MMGBSA Processing:')
        if '4' in flag:
            process4 = multiprocessing.Process(target=mmgbdec.main, args=(
                mpi_num, pdbcom_prmtop, pdbpro_prmtop, pdblig_prmtop, pdbcomsolvate_prmtop))
            process4.start()
            record.append(process4)
        if '5' in flag:
            process5 = multiprocessing.Process(target=nmode.main, args=(
                mpi_num, pdbcom_prmtop, pdbpro_prmtop, pdblig_prmtop, pdbcomsolvate_prmtop))
            process5.start()
            record.append(process5)

        elif flag == '6':  # 需要用户手动操作 无法多进程化
            dis_ang.main(traj, pdbcomsolvate_prmtop)
        else:
            raise RuntimeError('Invalid Code. Exit.')

        for process in record:
            process.join()  # 阻塞进程 等待任务结束

        print('\nAll Analysis Processing Complete.')
        print(''.center(80, '-'))


if __name__ == '__main__':
    pynalysis = Pynalysis()
    pynalysis.main()
