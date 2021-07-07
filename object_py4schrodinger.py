import csv
import getopt
import multiprocessing
import os
import re
import sys
import getopt
import time

from schrodinger import structure as struc
from schrodinger.application.glide import poseviewconvert as pvc
from schrodinger.job import jobcontrol as jc
from schrodinger.protein import getpdb

cwd = str(os.getcwd())
sys.path.append(cwd)
total = 0
now_complete = 0

def error_handler(error):
    print(error.__cause__)
    print(''.center(80, '-'), end='\n')


def success_handler(result):
    global total
    global now_complete
    now_complete += result
    print("%s Job(s) Successd / Total: %s" % (now_complete, total))


class Console:

    def __init__(self, pdbid='', ligname='') -> None:
        '''
        Script Core For Schrodinger Suite Analysis 
        Version 1.10

        Author: YH. W
        Last Update: 2021/07/06
        
        Property
        ----------
        pdbid: PDB ID
        ligname: 配体文件PATH
        '''

        self.pdbid = pdbid
        self.ligname = ligname
        self.pdbfile = ''           # PDB结构文件名
        self.minimized_file = ''    # Minimized化完成的文件名
        self.grid_file = ''         # 格点文件名
        self.lig_file = ''          # 内源配体文件名
        self.recep_file = ''        # 受体文件名
        self.dock_file = ''         # 对接结果文件名

    @staticmethod
    def load_st(st_file) -> object:
        '''
        读取结构

        Parameters
        ----------
        st_file: 需要读取结构的文件path 需要Schrodinger支持的文件格式

        Return
        ----------
        结构对象
        '''

        return next(struc.StructureReader(st_file))

    @staticmethod
    def __launch(cmd):
        '''
        使用jobcontrol启动一项job并等待结束

        Parameters
        ----------
        cmd： 等待执行的命令字符串
        '''

        cmd_list = cmd.split(' ')  # launch_job以列表形式提交参数
        job = jc.launch_job(cmd_list)
        print('JobId: %s' % job.JobId, end='\n')
        print('Job Name: %s' % job.Name, end='\n')
        print('Job Status: %s' % job.Status, end='\n')
        job.wait()  # 阻塞进程 等待Job结束

    @staticmethod
    def checkpdb(pdb:str):
        '''
        检查PDB ID合法性

        Return
        ----------
        合法返回True 否则返回False
        '''

        match = re.fullmatch(r'^\d[0-9a-zA-Z]{3,}$', pdb)
        if match:
            return True
        else:
            return False

    @staticmethod
    def check_ligname(ligname:str): #可能没有必要
        '''
        检查配体名称合法性

        Return
        ----------
        合法返回Match对象 否则返回False
        '''

        match = re.search('[0-9A-Z]{2,3}$', ligname)
        if match:
            return match
        else:
            return False

    def convert_format(self, file:str, suffix:str) -> str:
        '''
        mae/pdb格式转换为pdb/mae格式

        Parmeters
        ----------
        file: 需要转换的mae文件PATH
        suffix: 格式后缀名(pdb|mae)

        Reture
        ----------
        转换后的文件名
        '''

        st = self.load_st(file)
        if suffix == 'pdb':
            st.write(file.split('.')[0] + '.pdb')
            return file.split('.')[0] + '.pdb'
        elif suffix == 'mae':
            st.write(file.split('.')[0] + '.mae')
            return file.split('.')[0] + '.mae'

    def get_pdbid(self) -> str:
        '''
        获取用户输入的PDBID并检查合法性

        Return
        ----------
        PDB ID字符串

        '''

        if not self.pdbid:
            pdb = str(os.path.split(cwd)[1])
        else:
            pdb = str(self.pdbid)

        if self.checkpdb(pdb):
            self.pdbid = pdb    #[属性修改]修改PDB ID
            self.pdbfile = pdb + '.pdb' #[属性修改]修改原始结构文件PATH
            return pdb
        else:
            while True:
                pdb = str(
                    input('\n要自动获取PDB ID 请将晶体文件夹修改为PDB ID\n请手动输入晶体PDB ID:')).strip().upper()
                if self.checkpdb(pdb):
                    self.pdbid = pdb #[属性修改]修改PDB ID
                    self.pdbfile = pdb + '.pdb' #[属性修改]修改原始结构文件PATH
                    return pdb
                else:
                    print('请输入正确的PDB ID!')

    def get_ligname(self) -> str:
        '''
        尝试自动获取配体名 如有多个配体则获取用户输入的配体名并检查合法性

        Return
        ----------
        自动识别或手动输入的配体名称

        '''

        if self.pdbid:
            pdbid = self.pdbid
        else:
            pdbid = self.get_pdbid()

        if not os.path.exists(self.pdbfile):  # 下载PDB文件
            print('%s.pdb is not found. Would you like to download it? (Y/N)' % pdbid)
            __temp_flag = input().upper().strip()
            if __temp_flag == 'Y':
                getpdb.get_pdb(pdbid)
            else:
                print('No PDB File Found. Exit.')

        lis = os.popen(
            "cat %s.pdb | grep -w -E ^HET | awk '{print $2}'" % pdbid).readlines()   # 抓取pdb原始结构文件(不可是已处理过的结构)中单一entry关于小分子的描述
        lig = []

        for i in lis:
            passed = self.check_ligname(i.strip())  # 匹配得到配体候选列表
            if passed:
                lig_name = passed.group()
                if lig_name != 'HOH' and not lig_name in lig:  # 排除配体候选中的水分子与重复分子
                    lig.append(lig_name)

        if len(lig) == 1:
            ligname = str(lig[0])
            self.ligname = ligname  #[属性修改] 修改配体名称
            return ligname
        else:
            print('晶体%s包含多个配体:' % pdbid, ''.join(
                str(x)+' ' for x in lig), end='\n')
            while True:
                ligname = input('请手动指定配体名称:').strip().upper()
                if self.check_ligname(ligname) and ligname in lig:
                    self.ligname = ligname  #[属性修改] 修改配体名称
                    return ligname
                else:
                    print('请重新输入正确的配体名称!')
    

    def keep_chain(self, chain_name:str, pdbfile:str=None) -> str:
        '''
        读取PDB晶体文件并将单一链的结构输出为pdb文件

        Parameters
        ----------

        chain_name: 要保留的链名称
        pdbfile: 待处理原始文件

        Return
        ----------
        保留单链结构的文件名
        '''

        if not pdbfile:
            pdbfile = self.pdbfile

        st = self.load_st(pdbfile)  # 读取原始PDB结构
        st_chain_only = st.chain[chain_name].extractStructure()
        singlechain_file = '%s_chain_%s.mae' % (pdbfile.split('.')[0], chain_name)
        st_chain_only.write(singlechain_file)
        return singlechain_file

    def __preprocess(self) -> str:
        '''
        一些预处理操作：
        下载PDB文件 | 检查PDB晶体类型 Apo/单体/多聚体 | 是否保留单链   
        将会修改property：pdbfile

        Returen
        ----------
        处理完成的PDB结构文件名
        '''

        pdbid = self.get_pdbid()
        pdbfile = self.pdbfile

        if not os.path.exists(pdbfile):  # 下载PDB文件
            getpdb.get_pdb(pdbid)

        lig_lis = os.popen(
            "cat %s | grep -w -E ^HET | awk '{print $2}'" % pdbfile).readlines()

        if len(lig_lis) == 0:
            raise RuntimeError('%s为Apo蛋白晶体 无法自动处理.' % pdbid)

        elif len(lig_lis) > 1:
            print('\n')
            os.system('cat %s.pdb | grep -w -E ^HET' % pdbid)
            print('存在多个配体小分子 是否需要保留单链？(Y/N)')
            _flag = input().strip().upper() # 是否保留单链的标志

            if _flag == 'Y':
                chain = input('输入保留链名:').strip().upper()
                pdbfile = self.keep_chain(chain)
                self.pdbfile = pdbfile  # [属性修改]修改pdbfile为单链文件
                return pdbfile
            else:
                print('\nChoosed Intact Crystal.\n')
                return pdbfile  # 主要晶体文件仍为原始晶体文件
        else:
            return pdbfile
    
    def minimize(self, pdbfile:str=None) -> str:
        '''
        调用prepwizard模块自动优化PDB结构文件

        Parameters
        ----------
        pdbfile: 需要优化的文件

        Return
        ----------
        完成优化后的文件名

        '''

        if not pdbfile:
            pdbfile = self.__preprocess()  # 结构文件预处理

        minimized_file = str(pdbfile.split('.')[0] + '_minimized.mae')

        if os.path.exists(minimized_file):  # 如果已经进行过优化 为提高效率而跳过优化步骤
            self.minimized_file = minimized_file    # [属性修改] 修改最小化后的结构文件PATH
            return minimized_file

        prepwizard_command = 'prepwizard -f 3 -r 0.3 -propka_pH 7.0 -disulfides -s -j %s-Minimize %s %s' % (pdbfile.split('.')[0],
                                                                                                            pdbfile, minimized_file)
        self.__launch(prepwizard_command)   # 阻塞至任务结束

        if not os.path.exists(minimized_file):  # 判断Minimized任务是否完成(是否生成Minimized结束的结构文件)
            raise RuntimeError(
                '%s Crystal Minimization Process Failed' % pdbfile.split('.')[0])  # 无法被优化的晶体结构
        else:
            print('\nPDB Minimized File', minimized_file, 'Saved.\n')
            self.minimized_file = minimized_file    # [属性修改] 修改最小化后的结构文件PATH
            return minimized_file
    
    def __get_ligmol_info(self) -> str:
        '''
        以ligname为KEY 查找Maestro文件中的Molecule Number

        Return
        ----------
        Ligand 所在Molecule Number

        '''

        minimized_file = self.minimized_file
        ligname = self.ligname
        st = self.load_st(minimized_file)  # 载入结构对象

        def _get_mol(st):   # 获取结构中的配体所在Molecule对象
            residues = st.residue
            for res in residues:
                if res.pdbres.strip() == '%s' % ligname:
                    molnum = res.molecule_number
                    yield st.molecule[molnum]

        mol = next(_get_mol(st))

        if len(mol.residue) != 1:  # 判断该molecule是否仅包括小分子本身(是否存在共价连接)
            print('%s in %s 配体分子与残基可能存在共价连接 将尝试自动删除共价键' % (ligname, minimized_file))
            bonds = st.bond

            for bond in bonds:
                resname1 = bond.atom1.getResidue().pdbres.strip()
                resname2 = bond.atom2.getResidue().pdbres.strip()
                if resname1 == '%s' % ligname or resname2 == '%s' % ligname:
                    if resname1 != resname2:
                        bond_to_del = bond
            
            if not bond_to_del:
                raise RuntimeError('无法自动获取共价键')
            
            st.deleteBond(bond_to_del.atom1, bond_to_del.atom2)
            st.write(minimized_file)

            mol = next(_get_mol(st))

        return mol.number

    def grid_generate(self, pdbid:str=None, ligname:str=None, st_file:str=None, gridbox_size:int=20) -> str:
        '''
        自动编写glide grid输入文件并启动Glide Grid生成任务

        Parameters
        ----------
        pdbid: PDB ID
        ligname: 配体名称(RCSB ID)
        st_file: 需要构建grid box的文件PATH
        gridbox_size: grid box大小 默认20Å

        Return
        ----------
        生成的格点文件名

        '''
        
        lig_molnum = self.__get_ligmol_info()

        # 默认读取类属性中的信息
        if not pdbid:
            pdbid = self.pdbid
        if not ligname:
            ligname = self.ligname
        if not st_file:
            st_file = self.minimized_file
            
        outsize = gridbox_size + 10
        
        grid_file = pdbid + '_glide_grid_%s.zip' % ligname

        if os.path.exists(grid_file):   # 如果已经生成了格点文件则跳过生成过程
            self.grid_file = grid_file  # [属性修改] 修改格点文件PATH
            return grid_file

        with open('%s_grid_generate_%s.in' % (pdbid, ligname), 'w') as input_file:  # 编写glide输入文件
            input_file.write('GRIDFILE %s\n' % grid_file)  # 输出的Grid文件名
            input_file.write('INNERBOX 10,10,10\n')  # Box大小参数
            input_file.write('OUTERBOX %d,%d,%d \n' %
                            (outsize, outsize, outsize))  # Box大小参数
            input_file.write('LIGAND_MOLECULE %s\n' %
                            lig_molnum)  # 识别Ligand并设定grid box中心为质心
            input_file.write('RECEP_FILE %s' % st_file)  # 输入文件
        self.__launch('glide %s_grid_generate_%s.in -JOBNAME %s-%s-Grid-Generate' %
            (pdbid, ligname, pdbid, ligname))

        print('\nGrid File', grid_file, 'Saved.\n')
        self.grid_file = grid_file  # [属性修改] 修改格点文件PATH
        return grid_file

    def split_com(self, complex_file:str=None, ligname:str=None) -> list:
        '''
        拆分复合体为配体和受体

        Parameters
        ----------
        complex_file: 待拆分的复合物文件PATH
        ligname: 配体文件名(RCSB ID)

        Return
        ----------
        List: [ligfile_PATH, recepfile_PATH]

        '''

        if not ligname:
            if self.ligname:
                ligname = self.ligname
            else:
                ligname = self.get_ligname()
        if not complex_file:
            if self.minimized_file:
                complex_file = self.minimized_file  # 默认拆分Minimized完成的结构文件
            else:
                raise RuntimeError('No Complex File Provided.')

        pdbid = complex_file.split('_')[0]

        st = self.load_st(complex_file)
        residue_list = []
        for res in st.residue:
            if res.pdbres.strip() == '%s' % ligname:
                residue_list.append(res)

        comp = pvc.Complex(st, ligand_asl=residue_list[0].getAsl(), ligand_properties=st.property.keys())

        if len(residue_list) != 1:
            resname_list = []
            for res in residue_list:
                text = res.chain + ':' + res.pdbres.strip() + ' ' + str(res.resnum)
                resname_list.append(text)
            print('\nThere are %s "%s" in %s: ' % (len(residue_list), ligname, complex_file), resname_list, '\nThe First One is Selected')

        lig_file = '%slig_%s.mae' % (pdbid, ligname)
        self.lig_file = lig_file                 # [属性修改] 修改配体文件PATH
        recep_file = '%spro_%s.mae' % (pdbid, ligname)     
        self.recep_file = recep_file             # [属性修改] 修改受体文件PATH
        comp.writeLigand(lig_file)          # 生成并保存配体独立mae文件
        comp.writeReceptor(recep_file)      # 生成并保存受体独立mae文件

        st.write('%scom_%s.pdb' % (pdbid, ligname))
        self.convert_format(lig_file, 'pdb')        # 自动生成PDB格式
        self.convert_format(recep_file, 'pdb')

        return [lig_file, recep_file]

    def dock(self, pdbid:str=None, lig_file:str=None, grid_file:str=None, precision:str='SP', calc_rmsd:bool=False) -> str:
        '''
        一对一 多对一 glide dock任务输入文件编写与运行

        Parameters
        ----------
        pdbid: PDB ID
        lig_file: 配体文件PATH
        grid_file: 格点文件PATH
        precision: 对接精度(HTVS|SP|XP) 默认SP
        calc_rmsd: 是否计算rmsd to input ligand geometries 默认False
        '''

        if not pdbid:
            pdbid = self.pdbid
        if not lig_file:
            lig_file = self.lig_file
        if not grid_file:
            grid_file = self.grid_file

        print('Prepare to Docking...\n')
        lig_name = lig_file.split('.')[0].split('_')[-1]

        with open('%s_glide_dock_%s_%s.in' % (pdbid, lig_name, precision), 'w') as input_file:
            input_file.write('GRIDFILE %s\n' % grid_file)
            # 如果lig文件包含多个lig mol 将自动从一至末尾全部dock
            input_file.write('LIGANDFILE %s\n' % lig_file)
            input_file.write('PRECISION %s\n' % precision)  # HTVS SP XP
            if calc_rmsd == True:                           # 计算与原配体rmsd 默认False
                input_file.write('CALC_INPUT_RMS True\n')
            if precision == 'XP':
                input_file.write('WRITE_XP_DESC False\n')
                input_file.write('POSTDOCK_XP_DELE 0.5\n')

        self.__launch('glide %s_glide_dock_%s_%s.in -JOBNAME %s-Glide-Dock-%s-%s' %
            (pdbid, lig_name, precision, pdbid, lig_name, precision))

        c = os.system('mv %s-Glide-Dock-%s-%s_pv.maegz %s_glide_dock_%s_%s.maegz' %
                    (pdbid, lig_name, precision, pdbid, lig_name, precision))
        if c != 0:
            raise RuntimeError('%s-%s Gilde Docking Failed' % (pdbid, lig_name))

        print('\nDocking Result File:', '%s_glide_dock_%s_%s.maegz Saved.\n' %
            (pdbid, lig_name, precision))
        dock_file = '%s_glide_dock_%s_%s.maegz' % (pdbid, lig_name, precision)
        self.dock_file = dock_file                         #[属性修改] 修改对接结果文件PATH
        return dock_file

    def extra_data(self, path:str=None, ligname:str=None, precision:str='SP') -> dict:
        '''
        从对接完成的Maestro文件中提取数据

        Parameters
        ----------
        path: 对接完成的文件PATH
        ligname: 参与对接的配体名称
        precision: 已完成的对接工作精度

        Return
        ----------
        Dict
        Properties : Values

        '''

        file = struc.StructureReader(path)
        pdbid = self.pdbid
        if not path:
            path = self.dock_file
        if not ligname:
            ligname = self.ligname

        pro_st = next(file)
        lig_st = next(file)
        prop_dic = {}

        # 需要提取的Property
        prop_dic['PDB'] = pdbid
        prop_dic['Ligand'] = ligname
        prop_dic['Docking_Score'] = lig_st.property['r_i_docking_score']  # 对接分数
        prop_dic['rotatable_bonds'] = lig_st.property['i_i_glide_rotatable_bonds']
        prop_dic['ligand_efficiency'] = lig_st.property['r_i_glide_ligand_efficiency']
        prop_dic['evdw'] = lig_st.property['r_i_glide_evdw']
        prop_dic['ecoul'] = lig_st.property['r_i_glide_ecoul']
        prop_dic['energy'] = lig_st.property['r_i_glide_energy']
        prop_dic['einternal'] = lig_st.property['r_i_glide_einternal']
        prop_dic['emodel'] = lig_st.property['r_i_glide_emodel']
        prop_dic['precision'] = precision

        try:
            prop_dic['rmsd'] = lig_st.property['r_i_glide_rmsd_to_input']
        except KeyError:
            pass

        if precision == 'SP':
            prop_dic['lipo'] = lig_st.property['r_i_glide_lipo']
            prop_dic['hbond'] = lig_st.property['r_i_glide_hbond']
            prop_dic['metal'] = lig_st.property['r_i_glide_metal']
            prop_dic['rewards'] = lig_st.property['r_i_glide_rewards']
            prop_dic['erotb'] = lig_st.property['r_i_glide_erotb']
            prop_dic['esite'] = lig_st.property['r_i_glide_esite']

        elif precision == 'XP':
            prop_dic['XP_Hbond'] = lig_st.property['r_glide_XP_HBond']

        return prop_dic
        

class UI:

    @staticmethod
    def print_title():
        print(''.center(80,'*'),end='\n')
        print('Python Script For Schrodinger Suite Analysis'.center(80), end='\n')
        print(''.center(80,'*'),end='\n')

    def __init__(self, pdbid:str=None, ligname:str=None) -> None:
        '''
        User Interface Mode of Script
        '''
        console = Console()
        pdbid = console.get_pdbid()
        ligname = console.get_ligname()
        self.console = console

        self.print_title()
        print('\nPDB ID:', pdbid, end='\n')
        print('Entry Ligand:', ligname, end='\n')

    def __get_flag(self) -> str:

        while True:
            flag = input('''
    请输入需要进行的分析操作:

    1.PDB文件获取+优化
    2.PDB文件获取+优化+生成格点文件(Size 20Å)
    3.仅生成格点文件(自定义Size)
    4.自动执行PDB文件内源配体对接(SP精度)
    5.自动执行PDB文件内源配体对接(XP精度)
    6.非内源性配体对接

    0.退出

            ''')
            if re.match('^[0123456]$', flag):
                return flag
            else:
                print('请重新输入有效的操作代号!')
        
    def process(self) -> None:
        '''
        业务逻辑判断
        '''

        console = self.console
        flag = self.__get_flag()

        if flag == '0':
            sys.exit(0)

        elif len(flag) == 1 and flag in '12345':
            console.minimize()                  # 能量最小化
            split_lis = console.split_com()     # 拆分复合物
            lig_file = split_lis[0]             # 配体文件PATH获取

            if flag == '2':
                grid_file = console.grid_generate()     # 默认条件生成Grid文件
            
            elif flag == '3':
                size = int(input('Input Grid Box Size(Default 20Å):'))
                console.grid_generate(gridbox_size=size)    # 以指定Size生成Grid文件
            
            elif flag in '45':
                print('Ligand File:', lig_file)
                grid_file = console.grid_generate()
                print('Grid File:', grid_file)

                if flag == '4':
                    console.dock(precision='SP', calc_rmsd=True)
                else:
                    console.dock(precision='XP', calc_rmsd=True)

        elif flag == '6':
            lig_file = input('Input Ligand File PATH:').strip()
            precision = input('Input Docking Precision(HTVS|SP|XP):').strip().upper()

            lig_file = console.convert_format(lig_file, 'mae')
            console.minimize()
            grid_file = console.grid_generate()
            console.dock(lig_file=lig_file, precision=precision)

        else:
            raise RuntimeError('Wrong Input Code\nEXIT.')
            
        print(''.center(80, '-'), end='\n')

class Mutildock(Console):

    def __init__(self, pdbid, ligname) -> None:
        super().__init__(pdbid=pdbid, ligname=ligname)
        self.list_file = ''         # 单一基因的全部PDB晶体列表文件
        self.ligand_file = ''       # 外源配体文件
        self.precision = ''         # 对接精度
        self.cpus = ''              # 对接工作使用的进程数
        self.flag = ''              # 对接工作启动确认标记
        self.list_filename = ''     # 基因名称(晶体列表文件的名称)
        self.notpass = []           # 内源配体对接不成功的晶体列表
        self.dock_fail = []         # 外源配体对接不成功的晶体列表

    def usage():
        '''
Mutil-Dock Auto Processing.

Usage: 
UI Mode: run py4schrodinger.py
Mutil-dock Mode: run py4schrodinger.py -r|--receptor <receptors list file> [-l|--ligand <ligand file> -p|--precision <precision> -n|--cpu <cpus> -k|--no-check ]

Descrption
        Required
            -r, --receptor PATH of PDB ID List File(.txt) 
        Optional
            -l, --ligand PATH of Ligand File needed to be docked(.pdb)     
            -p, --precision Specify Docking Precision(HTVS|SP|XP)       
            -n, --cpu Specify Number of CPU used to dock     
            -k, --no-check Running Docking Without Check 'Y/N'      

        Note: Only crystals' own ligand but No exogenous ligand will be docked without [-l|--ligand].
'''

    def __process_argvs(self):
        '''
        解析命令行参数

        Parameters
        ----------
        argv: 命令行参数

        '''
        
        try:
            opts, argvs = getopt.getopt(sys.argv[1:], '-hkr:l:p:n:', ['help', 'no-check', 'receptor=', 'ligand=', 'precision=', 'cpu='])
        except getopt.GetoptError as err:
            print(str(err))
            print(self.usage.__doc__)
            sys.exit(1)

        for opt, arg in opts:       # 参数解析

            if opt in ('-h','--help'):
                print(self.usage.__doc__)
            elif opt in ('-r', '--receptor'):
                self.list_file = arg
                self.list_filename = self.list_file.split('.')[0].split(os.sep)[-1]
            elif opt in ('-l', '--ligand'):
                self.ligand_file = arg
            elif opt in ('-p', '--precision'):
                self.precision = arg
            elif opt in ('-k', '--no-check'):
                self.flag = 'Y'
            elif opt in ('-n', '--cpu'):
                self.cpus = int(arg)

    def __list_process(self) -> list:
        '''
        读取晶体列表文件信息

        '''
        pdb_list = []
        try:
            with open(self.list_file, 'r') as f:
                pdbs_withlig = f.readlines()        # 读取PDB列表中的每一行PDBID与配体名称信息
        except FileNotFoundError:
            print('Error: 未找到列表文件!')
            raise

        for i in pdbs_withlig:                      # 按逗号分割解析列表中的PDB ID与配体名称
            pdb = i.split(',')[0].strip().upper()   # PDB ID
            lig = i.split(',')[1].strip().upper()   # 配体名称
            pdb_list.append((pdb, lig))             # 将每一对作为元组储存至列表pdb_list
        
        return pdb_list

    def __check_before_dock(self, pdb_list:list) -> None:
        '''
        对接前检查
        '''

        for pdb, lig in pdb_list:  

            try:
                os.makedirs('./database/' + pdb)
            except FileExistsError:
                pass

            if not os.path.exists('./database/%s/%s.pdb' % (pdb, pdb)):  # 下载PDB文件
                getpdb.get_pdb(pdb)
                os.system('mv %s.pdb ./database/%s' % (pdb, pdb))

    def __check_fail(self, pdb_list: list) -> list:
        '''
        检查对接失败的晶体并从pdb_list中剔除

        Return
        ----------
        剔除了对接失败晶体的列表pdb_list
        '''

        notpass = []
        precision = self.precision

        for pdbid, ligname in pdb_list:  # 异常晶体跳过: APO & 共价键结合晶体 不会产生结果文件
            if not os.path.exists('./database/%s/%s_glide_dock_%s_%s.maegz' % (pdbid, pdbid, ligname, precision)):
                notpass.append((pdbid, ligname))
        if notpass:
            for not_exist, lig in notpass:
                pdb_list.remove((not_exist, lig))
            self.notpass = notpass
        
        return pdb_list


    def __autodock(self, pdbid:str, ligname:str, precision:str) -> int:
        '''
        自动化内源配体对接 for multidock

        Parameters
        ----------
        pdb: PDB ID字符串
        lig_name: 已核对的唯一配体
        precision: 对接精度

        '''
        os.chdir('./database/' + pdbid)

        cmd = '''cat %s.pdb | grep -w -E ^HET | awk '{if($2==\"%s\"){print $3}}' ''' % (pdbid, ligname)
        cmd_run = os.popen(cmd).readlines()

        if cmd_run:
            chain = re.match('[A-Z]', cmd_run[0].strip()).group()
        else:
            raise ValueError('No Match Ligand for %s' % ligname)

        if not chain:
            raise ValueError('No Chain Match')

        pdbfile = self.keep_chain(chain, '%s.pdb' % pdbid)

        print('\nPDB ID:', pdbid, end='\n')
        print('Entry Ligand:', ligname, end='\n')
        print('Chain: %s' % chain)

        # 如有Error将抛出异常 不再继续运行后续代码
        minimized_file = self.minimize(pdbfile)                 # 能量最小化

        lig_file = self.split_com(minimized_file, ligname)[0]   # 拆分复合物
        print('\nLigand File:', lig_file)

        grid_file = self.grid_generate(pdbid, ligname, minimized_file)    # 格点文件生成
        print('Grid File:', grid_file)

        self.dock(pdbid, lig_file, grid_file, precision, True)    # 对接 
        print('%s Self-Docking Job Complete.\n' % pdbid)
        print(''.center(80, '-'), end='\n')
        return 1                                                # 1项工作完成 返回以计数

    def __dock_multimode(self, pdbid:str, origin_ligname:str, ligand_file:str, precision:str) -> int:
        '''
        一对一 外源配体自动对接 for multidock

        Parameters
        ----------
        pdbid: 要对接到的受体PDB ID字符串
        origin_ligname: 对接位置原配体名称
        ligand_file: 要对接的外源配体文件PATH
        precision: 对接精度

        '''

        ligname = ligand_file.strip().split('.')[0]             # 文件名获取外源配体名
        ligand_file = self.convert_format(ligand_file, 'mae')   # 转换为Maestro格式

        os.chdir('./database/' + pdbid)                         # 切换工作目录
        ligand_file = '../../' + ligand_file                    # 重新定位外源配体文件PATH
        grid_file = '%s_glide_grid_%s.zip' % (pdbid, origin_ligname)

        print('\nPDB ID:', pdbid, end='\n')
        print('Ligand Name:', ligname)
        print('Grid File:', grid_file)

        print('Prepare to Docking...\n')

        # 撰写Dock输入文件
        with open('%s_glide_dock_%s_%s.in' % (ligname, origin_ligname, precision), 'w') as input_file:  
            input_file.write('GRIDFILE %s\n' % grid_file)       # 格点文件
            input_file.write('LIGANDFILE %s\n' % ligand_file)   # 外源配体PATH
            input_file.write('PRECISION %s\n' % precision)      # HTVS SP XP
            if precision == 'XP':
                input_file.write('WRITE_XP_DESC False\n')
                input_file.write('POSTDOCK_XP_DELE 0.5\n')

        # 与self.dock()不同点 需要区别同一位点不同配体的对接
        self.__launch('glide %s_glide_dock_%s_%s.in -JOBNAME %s-Glide-Dock-On-%s-%s-%s' %
            (ligname, origin_ligname, precision, ligname, pdbid, origin_ligname, precision))

        # 无论是否完成对接 返回码均为0 故通过是否产生对接结果文件进行对接成功与否区分        
        # 对接不成功时 无法产生对接结果文件
        if not os.path.exists('%s-Glide-Dock-On-%s-%s-%s_pv.maegz' % (ligname, pdbid, origin_ligname, precision)):
            raise RuntimeError('%s-%s-%s Gilde Dock Failed' % (ligname, pdbid, origin_ligname))

        # 重命名文件
        os.system('mv %s-Glide-Dock-On-%s-%s-%s_pv.maegz %s_glide_dock_on_%s_%s_%s.maegz' %
                (ligname, pdbid, origin_ligname, precision, ligname, pdbid, origin_ligname, precision))
        print('\nDocking Result File:', '%s_glide_dock_on_%s_%s_%s.maegz Saved.\n' %
            (ligname, pdbid, origin_ligname, precision))

        print('%s Docking on %s-%s Job Complete.\n' %(ligname, pdbid, origin_ligname))
        print(''.center(80, '-'), end='\n')
        return 1                                            # 1项对接工作完成 返回以计数

    def __auto_extradata(self, pdb_list:list) -> list:
        '''
        提取对接完成的数据构成字典 并将所有字典打包为列表返回

        Return
        ----------
        提取出的数据 list[dict1, dict2, ...]

        '''

        data = []                       # 总数据集
        dock_fail = []                  # 外源配体对接失败晶体列表
        precision = self.precision
        ligand_file = self.ligand_file

        for pdb, lig in pdb_list:

            origin_ligand = lig         # 原始配体名称
            prop_dic = self.extra_data('./database/%s/%s_glide_dock_%s_%s.maegz' %
                                    (pdb, pdb, origin_ligand, precision), origin_ligand)
            data.append(prop_dic)

            if ligand_file:             # 存在外源配体对接需求时
                ligname = ligand_file.strip().split('.')[0]
                dock_result_file = './database/%s/%s_glide_dock_on_%s_%s_%s.maegz' % (
                    pdb, ligname, pdb, origin_ligand, precision)
                if os.path.exists(dock_result_file):
                    ex_dic = self.extra_data(dock_result_file, ligname, precision)
                    data.append(ex_dic)
                else:                   # 对接不成功
                    dock_fail.append((pdb, origin_ligand, ligname))      

        self.dock_fail = dock_fail
        return data

    def __save_data(self, data):

        ligand_file = self.ligand_file
        list_filename = self.list_filename
        precision = self.precision
        notpass = self.notpass
        dock_fail = self.dock_fail

        # SP模式下与XP模式下产生对接结果项目不同
        prop_xp = ['PDB', 'Ligand', 'Docking_Score', 'rmsd', 'precision', 'ligand_efficiency', 'XP_Hbond', 'rotatable_bonds',
                    'ecoul', 'evdw', 'emodel', 'energy', 'einternal']
        prop_sp = ['PDB', 'Ligand', 'Docking_Score', 'rmsd', 'precision', 'ligand_efficiency', 'rotatable_bonds',
                    'ecoul', 'evdw', 'emodel', 'energy', 'einternal','lipo', 'hbond', 'metal', 'rewards', 'erotb', 'esite']
        withlig = '_' + ligand_file.split('.')[0]

        with open(list_filename + '_FINAL_RESULTS%s.csv' % withlig, 'w', encoding='UTF-8', newline='') as f:
            if precision == 'XP':
                writer = csv.DictWriter(f, fieldnames=prop_xp)
            elif precision == 'SP':
                writer = csv.DictWriter(f, fieldnames=prop_sp)
            writer.writeheader()    # 写入标头
            writer.writerows(data)  # 写入数据 自动匹配标头列
        
        if notpass:
            print('\nAbandoned Crystal(s): ', notpass)
            with open('./abandon.txt', 'a') as f:
                f.write('\n' + list_filename + '\n')
                for p, l in notpass:
                    f.write(p + ',' + l + '\n')
        
        if dock_fail:
            print('Docking Failed Crystal(s): ', dock_fail)

    def multidock(self) -> None:
        '''
        自动多进程处理多个PDB晶体并完成自动对接 提取对接结果数据并保存为CSV文件

        Parameters
        ----------
        pdb_list: （同一基因）包含多个pdb晶体ID的列表

        '''
        self.__process_argvs()              # 命令行参数解析
        pdb_list = self.__list_process()    # 列表文件解析
        self.__check_before_dock(pdb_list)  # 对接前检查
        UI.print_title()

        # 基本信息显示
        print('\nProcessing Input List...\n')
        print('All Crystals to be Processed:', ''.join(str(x)+' ' for x in pdb_list))
        global total
        total = len(pdb_list)
        print('\nTOTAL CRYSTALS:', total)

        if not self.flag:
            flag = input('\nContinue to Dock ? (Y/N)\n').strip().upper()
            self.flag = flag
        else:
            flag = self.flag

        if flag != 'Y':
            sys.exit(1)

        if not self.precision:
            precision = input('Input Precision of Docking(HTVS|SP|XP):').strip().upper()
            self.precision = precision
        else:
            precision = self.precision

        if not self.cpus:
            print('\nNumber of Total CPU:', multiprocessing.cpu_count(), end='\n')
            cpus = int(input('请输入需要使用的CPU核心数量:'))
            self.cpus = cpus

        print('Using Number of CPU: %s' % cpus)
        print('Docking Precision: %s' % precision)

        global now_complete
        now_complete = 0

        # 多进程对接
        multiprocessing.set_start_method('spawn')
        pool1 = multiprocessing.Pool(cpus, maxtasksperchild=1)  # 每个进程必须仅使用单线程
        pool2 = multiprocessing.Pool(cpus, maxtasksperchild=1)

        for pdbid, ligand in pdb_list:                          # 采用进程池控制多线程运行
            pool1.apply_async(self.__autodock, (pdbid, ligand, precision),
                            callback=success_handler, error_callback=error_handler)
            time.sleep(0.5)

        pool1.close()   # 进程池关闭 不再提交新任务
        pool1.join()    # 阻塞进程 等待全部子进程结束

        pdb_list = self.__check_fail(pdb_list)
        
        if self.ligand_file:                                    # 存在外源性配体对接需求
            ligand_file = self.ligand_file

            print('\n')
            print(''.center(80,'-'))
            print('User-Defined Ligand Docking'.center(80))
            print(''.center(80,'-'))
            print('\nUser-Defined Ligand:', ligand_file.split('.')[0])

            global now_complete
            now_complete = 0
            for pdb_code, lig in pdb_list:
                pool2.apply_async(self.__dock_multimode, (pdb_code, lig, ligand_file, precision,),
                                callback=success_handler, error_callback=error_handler)
                time.sleep(0.5)

            pool2.close()
            pool2.join()

        data = self.__auto_extradata(pdb_list)
        self.__save_data(data)

        print('\nAll Docking Jobs Done.\n')

def main():

    ui = UI()
    ui.process()


if __name__ == '__main__':
    main()