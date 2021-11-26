import os
import re
import sys

root_path = os.path.abspath(os.path.dirname(__file__)).split('src')[0]  # 项目路径 绝对路径
lib_path = root_path + 'lib' + os.sep                                   # 库文件夹路径
pdb_path = root_path.split('automatedMD')[0]                            # PDB项目绝对路径(如果有)
vsw_path = pdb_path + 'vsw' + os.sep                                    # VSW库文件夹路径
list_path = lib_path + 'list' + os.sep                                  # 存储各基因所含PDBID列表文件库路径
dockfiles_path = lib_path + 'dockfiles' + os.sep                        # 对接数据文件存储库路径

# Nuclear Receptor Gene列表 此后可更换为其他家族基因
gene_dict = {'NR1A1': 'TR_alpha',
             'NR1A2': 'TR_beta',
             'NR1B1': 'RAR_alpha',
             'NR1B2': 'RAR_beta',
             'NR1B3': 'RAR_gamma',
             'NR1C1': 'PPAR_alpha',
             'NR1C2': 'PPAR_beta',
             'NR1C3': 'PPAR_gamma',
             'NR1D1': 'Rev-erb',
             'NR1F1': 'ROR_alpha',
             'NR1F2': 'ROR_beta',
             'NR1F3': 'ROR_gamma',
             'NR1H3': 'LXR_alpha',
             'NR1H2': 'LXR_beta',
             'NR1H4': 'FXR_alpha',
             'NR1H5': 'FXR_beta',
             'NR1I1': 'VDR',
             'NR1I2': 'PXR',
             'NR1I3': 'CAR',
             'NR2A1': 'HNF4_alpha',
             'NR2A2': 'HNF4_gamma',
             'NR2B1': 'RXR_alpha',
             'NR2B2': 'RXR_beta',
             'NR2B3': 'RXR_gamma',
             'NR2C1': 'TR2',
             'NR2C2': 'TR4',
             'NR2E2': 'TLL',
             'NR2E3': 'PNR',
             'NR2F1': 'COUP-TFI',
             'NR2F2': 'COUP-TFII',
             'NR2F6': 'EAR2',
             'NR3A1': 'ER_alpha',
             'NR3A2': 'ER_beta',
             'NR3B1': 'ERR_alpha',
             'NR3B2': 'ERR_beta',
             'NR3B3': 'ERR_gamma',
             'NR3C1': 'GR',
             'NR3C2': 'MR',
             'NR3C3': 'PR',
             'NR3C4': 'AR',
             'NR4A1': 'NGFIB',
             'NR4A2': 'NURR1',
             'NR4A3': 'NOR1',
             'NR5A1': 'SF1',
             'NR5A2': 'LRH1',
             'NR6A1': 'GCNF',
             'NR0B1': 'DAX1',
             'NR0B2': 'SHP'}

gene_list = list(gene_dict.keys())

# 化合物库路径 应及时添加和更换
lib_natural1500 = '/home/omnisky/database/1500/natural1500-3d.maegz'
lib_bioactive4413 = '/home/omnisky/database/4413/4413-3D.mae'
lib_chemdiv50000 = '/home/omnisky/database/chemdiv50000/chemdiv50000.maegz'
lib_chembridge110000 = '/home/omnisky/database/chembridge110000/Chembridge110000.maegz'

try:
    from schrodinger.job import jobcontrol as jc
except ImportError:
    raise ImportError('No module named schrodinger.\nPlease run this script with "run py4schrodinger.py" or in the schrodinger.ve environment.')


class VSW:
        
    '''
Python Script for Virtual Screening Workflow.
    '''

    def __init__(self, gene='', jobname='') -> None:
        self.gene = gene                # 需要进行VSW的基因
        self.pdb_list = []              # 由基因索取到的PDB ID列表 由元组(PDBID, Ligand)组成
        self.grid_num = ''              # 总计需要进行VSW的PDB数量(创建的pipeline数量)
        self.jobname = jobname          # 任务名称

    @staticmethod
    def _launch(cmd):
        '''
        使用jobcontrol启动一项job并等待结束

        Parameters
        ----------
        cmd ： str
            等待执行的命令字符串
        '''

        cmd_list = cmd.split(' ')  # launch_job以列表形式提交参数
        job = jc.launch_job(cmd_list)
        print('JobId: %s' % job.JobId, end='\n')
        print('Job Name: %s' % job.Name, end='\n')
        print('Job Status: %s' % job.Status, end='\n')
        job.wait()  # 阻塞进程 等待Job结束

    def get_gene(self) -> str:
        '''
        通过用户输入获取需要进行VSW的基因

        Reture
        ----------
        str
            需要进行VSW的基因名称
        '''
        gene = input('Please Enter the Gene needed to VSW:').strip().upper()
        if gene in gene_list:
            self.gene = gene                                                    # 基因位于默认的NR基因列表中 可以自动获取列表
            return gene
        else:
            print('No %s Detected in Default List.\nPlease Provide PDB ID List Files.' % gene)

    def get_pdblist(self, gene:str=None) -> list:
        '''
        获取由基因指定的某一成员全部PDB ID(存储于 lib/list/ 下的.txt文件中)
        Parameter
        ----------
        gene : str
            指定的基因名称
        
        Return 
        ----------
        list
            由基因指定成员所有的 (PDBID, ligand)元组组成的列表
        '''
        if not gene:
            gene = self.gene

        pdb_list = []
        try:
            with open(list_path + '%s.txt' % gene) as f:
                _temp_list = f.read().splitlines()                                  # 读取文件内容
                for _temp_text in _temp_list:
                    pdb = _temp_text.split(',')[0]                                  # PDB ID    
                    lig = _temp_text.split(',')[1]                                  # 配体名称
                    pdb_list.append((pdb, lig))
        except FileNotFoundError:
            raise FileNotFoundError('No list file of %s found in list directoy.' % gene)

        print('\nGene: %s' % gene)
        print('Gene Abbreviation: %s' % gene_dict[gene])
        print('PDB(s) of Gene %s: %s' % (gene, pdb_list))

        self.pdb_list = pdb_list
        return pdb_list

    def _file_check(self, pdb_list:list=None):
        '''
        检查所需文件是否存在

        Parameter
        ----------
        pdb_list : list
            需要进行检查的PDB ID列表
        
        '''
        if not pdb_list:
            pdb_list = self.pdb_list
            
        for pdb, lig in pdb_list:
            grid_file = '%s_glide_grid_%s.zip' % (pdb, lig)
            if not os.path.exists(dockfiles_path + pdb + os.sep + grid_file):
                raise FileNotFoundError('%s is not found. Please Generate Grid File Firstly.' % grid_file)
    
    def get_grid_path(self, pdb_list:list=None) -> list:
        '''
        获取所有PDB的格点文件路径
        
        Parameter
        ---------
        pdb_list : list
            PDB列表

        Return
        ---------
        list
            包含所有PDB的格点文件绝对路径的列表
        '''
        if not pdb_list:
            pdb_list = self.pdb_list

        grid_path_list = []

        for pdb, lig in pdb_list:
            grid_file_path = dockfiles_path + pdb + os.sep + '%s_glide_grid_%s.zip' % (pdb, lig)
            grid_path_list.append(os.path.abspath(grid_file_path))
        
        return grid_path_list

    def _get_lig_file(self):
        '''
        确定用于VSW的化合物库文件

        Return
        ----------
        str
            化合物库文件
        '''
        while True:
            _ligand_flag = input('''
Enter the code of Ligand Database to be used:
    
    1.  Natural Compounds 1500 Database
    2.  Bioactivate Compounds 4413 Database
    3.  Chemdiv 50000 Compounds Database
    4.  Chembridge 110000 Compounds Database

    5.  Specified Other Database
    
    0.  Exit
    ''')
            if not re.match('^[012345]$', _ligand_flag):
                print('Invalid Code, Please Try Again.')
            else:
                break
        
        if _ligand_flag == '1':
            lig_file = lib_natural1500
            lig_name = '1500'
        elif _ligand_flag == '2':
            lig_file = lib_bioactive4413
            lig_name = '4413'
        elif _ligand_flag == '3':
            lig_file = lib_chemdiv50000
            lig_name = '50000'
        elif _ligand_flag == '4':
            lig_file = lib_chembridge110000
            lig_name = '110000'
        elif _ligand_flag == '5':
            lig_file = input('Enter the PATH of Database File(.mae|.maegz):').strip()
            lig_name = os.path.basename(lig_file).split('.')[0]
        else:
            print('Invalid Code, Exit.')
            sys.exit(1)

        if not os.path.exists(lig_file):
            raise FileNotFoundError('No such file or directory: %s' % lig_file)

        print('\nSelected Database File: %s\n' % lig_file)
        self.lig_file = lig_file
        self.lig_name = lig_name

        if not self.jobname:
            self.jobname = self.gene + '_' + self.lig_name + '_VSW'

        return lig_file

    def gen_input_file(self, pdb_list:list, lig_file:str):
        '''
        生成vsw输入文件

        Parameter
        ---------
        pdb_list : list
            将要作为受体进行VSW的所有PDB ID
        lig_file : str
            用于VSW的化合物库文件PATH

        Return
        ----------
        str
            生成的inp输入文件名
        '''
        self._file_check(pdb_list)                                                  # 检查文件存在
        grid_path_list = self.get_grid_path(pdb_list)                                    # 生成Grid文件名列表

        if not self.jobname:
            jobname = 'VSW'
        else:
            jobname = self.jobname
        
        self.grid_num = len(grid_path_list)                                               # 即将为输入文件写入的Pipeline数量 即PDB总数

        input_file = jobname + '.inp'
        with open(vsw_path + input_file, 'w') as f:                                      # VSW 输入文件编写
            f.write('# Run as: $SCHRODINGER/vsw <inputfile> \n')

            # 载入配体数据集
            f.write('\n[SET:ORIGINAL_LIGANDS]\n')                                     
            f.write('    VARCLASS   Structures\n')
            f.write('    FILES   %s,\n' % lig_file)                                    # 化合物库PATH

            # 载入所有受体数据集
            for grid_path in grid_path_list:
                index = grid_path_list.index(grid_path) + 1                         # 默认从0开始索引
                f.write('\n[SET:GRID_%s]\n' % index)
                f.write('    VARCLASS   Grid\n')
                f.write('    FILE   %s\n' % grid_path)
            
            # 执行步骤
            for n in range(1, self.grid_num + 1):

                # PRE_DOCK_HTVS 参数                                   
                f.write('\n[STAGE:PRE_DOCK_HTVS_%s]\n' % n)                         
                f.write('    STAGECLASS   gencodes.RecombineStage\n')
                f.write('    INPUTS   ORIGINAL_LIGANDS,\n')
                f.write('    OUTPUTS   RECOMBINE_OUT_%s,\n' % n)
                f.write('    NUMOUT   njobs\n')
                f.write('    OUTFORMAT   maegz\n')
                f.write('    MIN_SUBJOB_STS   4000\n')
                f.write('    MAX_SUBJOB_STS   40000\n')
                f.write('    GENCODES   YES\n')
                f.write('    OUTCOMPOUNDFIELD   s_vsw_compound_code\n')
                f.write('    OUTVARIANTFIELD   s_vsw_variant\n')
                f.write('    UNIQUEFIELD   NONE\n')

                # DOCK_HTVS 参数
                f.write('\n[STAGE:DOCK_HTVS_%s]\n' % n)
                f.write('    STAGECLASS   glide.DockingStage\n')
                f.write('    INPUTS   RECOMBINE_OUT_%s, GRID_%s\n' % (n, n))
                f.write('    OUTPUTS   HTVS_OUT_%s,\n' % n)
                f.write('    RECOMBINE   NO\n')
                f.write('    PRECISION   HTVS\n')
                f.write('    UNIQUEFIELD   s_vsw_compound_code\n')
                f.write('    PERCENT_TO_KEEP   10.0\n')                 # 保留前10%
                f.write('    DOCKING_METHOD   confgen\n')
                f.write('    POSES_PER_LIG   1\n')
                f.write('    BEST_BY_TITLE   NO\n')
                f.write('    LIG_VSCALE   0.8\n')
                f.write('    LIG_CCUT   0.15\n')
                f.write('    MAXATOMS   300\n')
                f.write('    MAXROTBONDS   50\n')
                f.write('    AMIDE_MODE   penal\n')
                f.write('    POSE_OUTTYPE   LIB\n')
                f.write('    POSTDOCK   YES\n')
                f.write('    POSTDOCKSTRAIN   NO\n')
                f.write('    COMPRESS_POSES   YES\n')
                f.write('    EPIK_PENALTIES   YES\n')
                f.write('    FORCEPLANAR   NO\n')

                # PULL_HTVS 参数
                f.write('\n[STAGE:PULL_HTVS_%s]\n' % n)
                f.write('    STAGECLASS   pull.PullStage\n')
                f.write('    INPUTS   HTVS_OUT_%s, RECOMBINE_OUT_%s\n' % (n, n))
                f.write('    OUTPUTS   HTVS_OUT_ORIG_%s,\n' % n)
                f.write('    UNIQUEFIELD   s_vsw_compound_code\n')

                # PRE_DOCK_SP 参数
                f.write('\n[STAGE:PRE_DOCK_SP_%s]\n' % n)
                f.write('    STAGECLASS   gencodes.RecombineStage\n')
                f.write('    INPUTS   HTVS_OUT_ORIG_%s,\n' % n)
                f.write('    OUTPUTS   DOCK_SP_%s_INPUT,\n' % n)
                f.write('    NUMOUT   njobs\n')
                f.write('    OUTFORMAT   maegz\n')
                f.write('    MIN_SUBJOB_STS   300\n')
                f.write('    MAX_SUBJOB_STS   3000\n')
                f.write('    GENCODES   NO\n')
                f.write('    UNIQUEFIELD   s_vsw_compound_code\n')

                # DOCK_SP 参数
                f.write('\n[STAGE:DOCK_SP_%s]\n' % n)
                f.write('    STAGECLASS   glide.DockingStage\n')
                f.write('    INPUTS   DOCK_SP_%s_INPUT, GRID_%s\n' % (n, n))
                f.write('    OUTPUTS   SP_OUT_%s,\n' % n)
                f.write('    RECOMBINE   NO\n')
                f.write('    PRECISION   SP\n')
                f.write('    UNIQUEFIELD   s_vsw_compound_code\n')
                f.write('    PERCENT_TO_KEEP   50\n')                              # 保留前50%
                f.write('    DOCKING_METHOD   confgen\n')
                f.write('    POSES_PER_LIG   1\n')
                f.write('    WRITE_XP_DESC   NO\n')
                f.write('    NENHANCED_SAMPLING   1\n')
                f.write('    BEST_BY_TITLE   NO\n')
                f.write('    LIG_VSCALE   0.8\n')
                f.write('    LIG_CCUT   0.15\n')
                f.write('    MAXATOMS   300\n')
                f.write('    MAXROTBONDS   50\n')
                f.write('    AMIDE_MODE   penal\n')
                f.write('    POSE_OUTTYPE   LIB\n')
                f.write('    POSTDOCK   YES\n')
                f.write('    POSTDOCKSTRAIN   NO\n')
                f.write('    COMPRESS_POSES   YES\n')
                f.write('    EPIK_PENALTIES   YES\n')
                f.write('    FORCEPLANAR   NO\n')

                # PULL_SP 参数
                f.write('\n[STAGE:PULL_SP_%s]\n' % n)
                f.write('    STAGECLASS   pull.PullStage\n')
                f.write('    INPUTS   SP_OUT_%s, HTVS_OUT_ORIG_%s\n' % (n, n))
                f.write('    OUTPUTS   SP_OUT_ORIG_%s,\n' % n)
                f.write('    UNIQUEFIELD   s_vsw_variant\n')

                # PRE_DOCK_XP 参数
                f.write('\n[STAGE:PRE_DOCK_XP_%s]\n' % n)
                f.write('    STAGECLASS   gencodes.RecombineStage\n')
                f.write('    INPUTS   SP_OUT_ORIG_%s,\n' % n)
                f.write('    OUTPUTS   DOCK_XP_%s_INPUT,\n' % n)
                f.write('    NUMOUT   njobs\n')
                f.write('    OUTFORMAT   maegz\n')
                f.write('    MIN_SUBJOB_STS   20\n')
                f.write('    MAX_SUBJOB_STS   200\n')
                f.write('    GENCODES   NO\n')
                f.write('    UNIQUEFIELD   s_vsw_compound_code\n')

                # DOCK_XP 参数
                f.write('\n[STAGE:DOCK_XP_%s]\n' % n)
                f.write('    STAGECLASS   glide.DockingStage\n')
                f.write('    INPUTS   DOCK_XP_%s_INPUT, GRID_%s\n' % (n, n))
                f.write('    OUTPUTS   XP_OUT_%s,\n' % n)
                f.write('    RECOMBINE   NO\n')
                f.write('    PRECISION   XP\n')
                f.write('    UNIQUEFIELD   s_vsw_compound_code\n')
                f.write('    PERCENT_TO_KEEP   20.0\n')                           # 保留前20%
                f.write('    DOCKING_METHOD   confgen\n')
                f.write('    POSES_PER_LIG   1\n')
                f.write('    WRITE_XP_DESC   NO\n')
                f.write('    BEST_BY_TITLE   YES\n')
                f.write('    LIG_VSCALE   0.8\n')
                f.write('    LIG_CCUT   0.15\n')
                f.write('    MAXATOMS   300\n')
                f.write('    MAXROTBONDS   50\n')
                f.write('    AMIDE_MODE   penal\n')
                f.write('    POSE_OUTTYPE   PV\n') 
                f.write('    POSTDOCK   YES\n')
                f.write('    POSTDOCKSTRAIN   NO\n')
                f.write('    COMPRESS_POSES   YES\n')
                f.write('    EPIK_PENALTIES   YES\n')
                f.write('    FORCEPLANAR   NO\n')

                # MMGBSA 参数
                f.write('\n[STAGE:MMGBSA_%s]\n' % n)
                f.write('    STAGECLASS   prime.MMGBSAStage\n')
                f.write('    INPUTS   XP_OUT_%s,\n' % n)
                f.write('    OUTPUTS   MMGBSA_%s,\n' % n)
        
            # 数据合并
            f.write('\n[STAGE:DOCKMERGE]\n')
            f.write('    STAGECLASS   glide.MergeStage\n')
            f.write('    INPUTS   ')
            for n in range(1, self.grid_num + 1):
                f.write('MMGBSA_%s,' % n)
            f.write('\n')
            f.write('    OUTPUTS   FINAL_DOCK_OUT,\n')
            f.write('    MAXPERLIG   1\n')
            f.write('    NREPORT   10000\n')
            f.write('    UNIQUEFIELD   s_vsw_compound_code\n')
            f.write('    OFFSETS   ')
            for n in range(1, self.grid_num + 1):
                f.write('0.0, ')
            f.write('\n')
            f.write('    MARKERS   ')
            for n in range(1, self.grid_num + 1):
                f.write('%s, ' % n)
            f.write('\n')

            # 数据输出
            f.write('\n[USEROUTS]\n')
            f.write('    USEROUTS   ')
            for n in range(1, self.grid_num + 1):
                f.write('XP_OUT_%s, MMGBSA_%s, ' % (n, n))
            f.write('FINAL_DOCK_OUT\n')
            f.write('    STRUCTOUT   FINAL_DOCK_OUT')

        return input_file
        
    def run(self):
        '''
        启动VSW
        '''
        if not self.jobname:
            jobname = 'VSW'
        else:
            jobname = self.jobname
        input_file = jobname + '.inp'
        input_file_path = os.path.abspath(vsw_path + input_file)
        if not os.path.exists(input_file_path):   # 检查文件是否存在
            raise FileNotFoundError('Input file %s is not found' % input_file)

        self._launch('vsw %s -OVERWRITE -host_glide localhost:48 -host_prime localhost:48 -adjust -NJOBS 48 -TMPLAUNCHDIR' % input_file)
        print('\nVSW job completed.\n')
        print(''.center(80,'-'))

class UI:
    '''
    用户交互界面
    '''

    @staticmethod
    def print_title():
        '''
        显示程序标题头
        '''
        print(''.center(80,'*'),end='\n')
        print('Python Script For Virtual Screening Workflow'.center(80), end='\n')
        print(''.center(80,'*'),end='\n')

    def __init__(self) -> None:

        vsw = VSW()
        try:
            os.chdir(vsw_path)
        except FileNotFoundError:
            os.makedirs(vsw_path)
            os.chdir(vsw_path)

        self.vsw = vsw
        self.print_title()

    def _get_flag(self) -> str:

        while True:
            flag = input('''
Please enter the code of analysis to be performed:

    1.  Virtual Screening in Nuclear Receptor Superfamily Default Database
    2.  Virtual Screening in Specified Database

    0.  Exit

            ''')
            if re.match('^[012]$', flag):
                return flag
            else:
                print('Invalid Code, Please Try Again.')
        
    def process(self) -> None:
        '''
        业务逻辑判断
        '''
        vsw = self.vsw
        _flag = self._get_flag()

        if _flag == '0':
            sys.exit(0)
        
        elif _flag == '1':
            gene = vsw.get_gene()                                   # 获取基因
            lig_file = vsw._get_lig_file()                          # 获取化合物库文件路径
            pdb_list = vsw.get_pdblist()                            # 获取PDB列表
            vsw._file_check(pdb_list)                               # 检查必要文件
            grid_path_list = vsw.get_grid_path(pdb_list)            # 获取格点文件路径

            _start_flag = input('Start VSW work?(Y/N):').strip().upper()
            if _start_flag == 'Y':
                input_file = vsw.gen_input_file(pdb_list, lig_file)     # 生成输入文件
                vsw.run()                                               # 开始VSW作业
        
        print(''.center(80, '-'), end='\n')

if __name__ == '__main__':
    if len(sys.argv) == 1:  # UI模式
        ui = UI()
        ui.process()
        







    

        


