import csv
import getopt
import multiprocessing
import os
import re
import sys
import getopt

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

    def __init__(self, pdbid='', ligname='', flag='') -> None:
        '''
        Property
        ----------
        pdbid: PDB ID
        ligname: 配体文件PATH
        flag: 处理模式
        '''

        self.pdbid = pdbid
        self.ligname = ligname
        self.flag = flag
        self.pdbfile = ''           # PDB结构文件PATH
        self.minimized_file = ''    # Minimized化完成的文件PATH
        self.grid_file = ''         # 格点文件PATH
        self.lig_file = ''          # 配体文件PATH
        self.recep_file = ''        # 受体文件PATH
        self.dock_file = ''         # 对接结果文件PATH

    @staticmethod
    def load_st(st_file):
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
    def checkpdb(pdb):
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
    def check_ligname(ligname): #可能没有必要
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

    def convert_format(self, file, suffix):
        '''
        mae/pdb格式转换为pdb/mae格式

        Parmeters
        ----------
        file: 需要转换的mae文件PATH
        suffix: 格式后缀名(pdb|mae)

        '''

        st = self.load_st(file)
        if suffix == 'pdb':
            st.write(file.split('.')[0] + '.pdb')
        elif suffix == 'mae':
            st.write(file.split('.')[0] + '.mae')

    def get_pdbid(self):
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

    def get_ligname(self):
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
    

    def keep_chain(self, chain_name, pdbfile=''):
        '''
        读取PDB晶体文件并将单一链的结构输出为pdb文件

        Parameters
        ----------

        chain_name: 要保留的链名称
        pdbfile: 待处理原始文件

        Return
        ----------
        保留单链结构的文件PATH
        '''

        if not pdbfile:
            pdbfile = self.pdbfile

        st = self.load_st(pdbfile)  # 读取原始PDB结构
        st_chain_only = st.chain[chain_name].extractStructure()
        singlechain_file = '%s_chain_%s.mae' % (pdbfile.split('.')[0], chain_name)
        st_chain_only.write(singlechain_file)
        return singlechain_file

    def __preprocess(self):
        '''
        一些预处理操作：
        下载PDB文件 | 检查PDB晶体类型 Apo/单体/多聚体 | 是否保留单链

        Returen
        ----------
        处理完成的PDB结构文件PATH
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
    
    def minimize(self, pdbfile=None):
        '''
        调用prepwizard模块自动优化PDB结构文件

        Parameters
        ----------
        pdb_file: 需要优化的文件

        Return
        ----------
        完成优化后的文件PATH

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
    
    def __get_ligmol_info(self):
        '''
        以ligname为KEY 查找Maestro文件中的Molecule Number

        Parameters
        ----------
        minimized_file: maestro可读文件PATH
        lig_name: 配体文件名称

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

    def grid_generate(self, pdbid=None, ligname=None, st_file=None, gridbox_size=20):
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

    def split_com(self, complex_file=None, ligname=None):
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

    def dock(self, pdbid=None, lig_file=None, grid_file=None, precision='SP', calc_rmsd=False):
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

    def extra_data(self, path=None, ligname=None, precision='SP'):
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
        

class ui:

    def get_flag():

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

    def print_title():
        print('\n')
        print('Python Script For Schrodinger Suite Analysis'.center(80, '-'), end='\n')



def main():
    console = Console()
    console.get_pdbid()
    console.get_ligname()
    console.minimize()
    console.grid_generate()
    console.split_com()
    console.dock()
