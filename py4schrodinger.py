# Please Run this Script with command 'run py4schrodinger.py'
#!/home/omnisky/schrodinger2020-3/run
'''

Python Script For Schrodinger Suite Analysis 
Version 1.04

Author YH. W
Last Update: 2021/05/19

'''

from getopt import GetoptError, getopt
import os
import re
import sys
import multiprocessing
import csv
import getopt
import time
from schrodinger.protein import getpdb
from schrodinger.job import jobcontrol as jc
from schrodinger import structure as struc
from schrodinger.application.glide import poseviewconvert as pvc
cwd = str(os.getcwd())
sys.path.append(cwd)


def error_handler(error):
    print(error.__cause__)


def launch(cmd):
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
    print('Job Running...')
    job.wait()  # 阻塞进程 等待Job结束
    print('\nJob %s Complete.\n' % job.Name)


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


def convert_mae_to_pdb(file):
    '''
    mae格式转换为pdb格式

    Parmeters
    ----------
    file: 需要转换的mae文件PATH

    '''
    st = load_st(file)
    st.write(file.split('.')[0]+'.pdb')


def convert_pdb_to_mae(file):
    '''
    pdb格式转换为mae格式

    Parmeters
    ----------
    file: 需要转换的pdb文件PATH

    '''
    st = load_st(file)
    st.write(file.split('.')[0]+'.mae')


def check_ligname(ligname):
    '''
    检查配体名称合法性
    '''

    return re.search('[0-9A-Z]{2,3}$', ligname)


def checkpdb(pdb):
    '''
    检查PDB ID合法性
    '''
    return re.fullmatch(r'^\d[0-9a-zA-Z]{3,}$', pdb)


def get_pdbid():
    '''
    获取用户输入的PDBID并检查合法性

    Return
    ----------
    PDB ID字符串

    '''

    pdb = str(os.path.split(cwd)[1])
    if checkpdb(pdb):
        return pdb
    else:
        while True:
            pdb = str(
                input('\n要自动获取PDB ID 请将晶体文件夹修改为PDB ID\n请手动输入晶体PDB ID:')).strip().upper()
            if checkpdb(pdb):
                return pdb
            else:
                print('请输入正确的PDB ID!')


def get_ligname(pdb_code):
    '''
    尝试自动获取配体名 如有多个配体则获取用户输入的配体名并检查合法性

    Parameters
    ----------
    pdb_code: PDB ID字符串

    Return
    ----------
    自动识别或手动输入的配体名称

    '''

    # 抓取pdb文件中单一entry关于小分子的描述
    lis = os.popen(
        "cat %s.pdb | grep -w -E ^HET | awk '{print $2}'" % pdb_code).readlines()
    lig = []

    for i in lis:
        passed = check_ligname(i.strip())  # 匹配得到配体候选列表
        if passed:
            lig_name = passed.group()
            if lig_name != 'HOH' and not lig_name in lig:  # 排除配体候选中的水分子与重复分子
                lig.append(lig_name)

    if len(lig) == 1:
        return lig[0]
    else:
        print('晶体%s包含多个配体:' % pdb_code, ''.join(
            str(x)+' ' for x in lig), end='\n')
        while True:
            ligname = input('请手动指定配体名称:').strip().upper()
            if check_ligname(ligname) and ligname in lig:
                return ligname
            else:
                print('请重新输入正确的配体名称!')


def keep_chain(pdb_code, origin_file, chain_name):
    '''
    读取输入结构并将单一链的结构输出为pdb文件

    Parameters
    ----------
    pdb_code: PDB ID字符串
    origin_file: 待处理原始文件
    chain_name: 要保留的链名称

    Return
    ----------
    保留单链结构的文件PATH
    '''

    st = load_st(origin_file)  # 读取原始结构
    st_chain_only = st.chain[chain_name].extractStructure()
    file = '%s_chain_%s.mae' % (pdb_code, chain_name)
    st_chain_only.write(file)
    return file


def minimized(pdb_code, pdb_file):
    '''
    调用prepwizard模块自动优化

    Parameters
    ----------
    pdb_code: PDB ID字符串
    pdb_file: 需要优化的文件

    Return
    ----------
    完成优化后的文件PATH

    '''

    minimized_file = pdb_code + '_minimized.mae'

    if os.path.exists(minimized_file):
        return minimized_file

    prepwizard_command = 'prepwizard -f 3 -r 0.3 -propka_pH 7.0 -disulfides -s -j %s-Minimize %s %s' % (pdb_code,
                                                                                                        pdb_file, minimized_file)
    launch(prepwizard_command)
    print('\nPDB Minimized File', minimized_file, 'Saved.\n')
    return minimized_file


def get_lig_info(minimized_file, lig_name):
    '''
    以lig_name为KEY 查找Maestro文件中的分子&链&残基 对象
    并打包为元组返回

    Parameters
    ----------
    minimized_file: maestro可读文件PATH
    lig_name: 配体文件名称

    Return
    ----------
    tuple: (molecule, chain, residue)

    '''

    cmd = os.popen("cat %s | grep %s | awk '{print $1}'" % (
        minimized_file, lig_name))

    lig_atomnum = cmd.readlines()[2].strip()  # 如果有多个同名配体小分子 默认抓取到首个小分子及其所在链

    # 从结构中自动抓取ligand所在的molecule number(非residue number)
    st = load_st(minimized_file)  # 载入结构对象
    atom = st.atom[int(lig_atomnum)]  # 配体包含的原子对象
    mol = atom.getMolecule()  # 配体分子对象
    chain = atom.getChain()  # 配体所在的链对象
    residue = atom.getResidue()

    if len(mol.residue) != 1:  # 判断该molecule是否仅包括小分子本身(是否存在共价连接)
        raise RuntimeError('%s in %s 配体分子与蛋白残基可能存在共价连接 请手动删除共价键后重试\n' % (lig_name, minimized_file))

    return (mol, chain, residue)


def grid_generate(pdb_code, lig_name, minimized_file, gridbox_size=20):
    '''
    自动编写glide grid输入文件并启动Glide Grid生成任务

    Parameters
    ----------
    pdb_code: PDB ID
    lig_name: 配体名称字符串
    minimized_file: 需要构建grid box的文件PATH
    gridbox_size: grid box大小 默认20Å

    Return
    ----------
    生成的格点文件名

    '''

    lig_info = get_lig_info(minimized_file, lig_name)
    lig_molnum = lig_info[0].number

    size = gridbox_size + 10
    grid_file = pdb_code + '_glide_grid.zip'

    with open('%s_grid_generate.in' % pdb_code, 'w') as input_file:  # 编写glide输入文件
        input_file.write('GRIDFILE %s\n' % grid_file)  # 输出的Grid文件名
        input_file.write('INNERBOX 10,10,10\n')  # Box大小参数
        input_file.write('OUTERBOX %d,%d,%d \n' %
                         (size, size, size))  # Box大小参数
        input_file.write('LIGAND_MOLECULE %s\n' %
                         lig_molnum)  # 识别Ligand并设定grid box中心为质心
        input_file.write('RECEP_FILE %s_minimized.mae' % pdb_code)  # 输入文件
    launch('glide %s_grid_generate.in -JOBNAME %s-Grid-Generate' %
           (pdb_code, pdb_code))

    print('\nGrid File', grid_file, 'Saved.\n')
    return grid_file


def split_com(pdb_code, lig_name, complex_file):
    '''
    拆分复合体为配体和受体

    Parameters
    ----------
    pdb_code: PDB ID字符串
    lig_name: 配体名
    complex_file: 待拆分的复合物文件PATH

    Return
    ----------
    包含配体与受体文件名的列表

    '''

    st = load_st(complex_file)

    try:
        comp = pvc.Complex(
            st, ligand_asl='res. %s' % lig_name, ligand_properties=st.property.keys()  # 同一条链只有一个配体分子
        )
    except RuntimeError:
        resnum = os.popen("cat %s_minimized.mae | grep %s | awk '{print $6}'" % (  # 同一条链有多个同名配体分子
            pdb_code, lig_name)).readlines()[2].strip()
        comp = pvc.Complex(st, ligand_asl='res.num %s' %
                           resnum, ligand_properties=st.property.keys())
    lig_file = '%slig.mae' % pdb_code
    recep_file = '%spro.mae' % pdb_code
    comp.writeLigand(lig_file)  # 生成并保存配体独立mae文件
    comp.writeReceptor(recep_file)  # 生成并保存受体独立mae文件

    st.write('%scom.pdb' % pdb_code)
    convert_mae_to_pdb(lig_file)  # 自动生成PDB格式
    convert_mae_to_pdb(recep_file)

    return [lig_file, recep_file]


def dock(pdb_code, lig_file, grid_file, precision='SP', calc_rmsd=False):
    '''
    一对一 多对一 glide dock任务输入文件编写与运行

    Parameters
    ----------
    pdb_code: PDB ID字符串
    lig_file: 配体文件PATH
    grid_file: 格点文件PATH
    precision: 对接精度(HTVS|SP|XP) 默认SP
    calc_rmsd: 是否计算rmsd to input ligand geometries 默认False

    '''
    print('Prepare to Docking...\n')
    with open('%s_glide_dock_%s.in' % (pdb_code, precision), 'w') as input_file:
        input_file.write('GRIDFILE %s\n' % grid_file)
        # 如果lig文件包含多个lig mol 将自动从一至末尾全部dock
        input_file.write('LIGANDFILE %s\n' % lig_file)
        input_file.write('PRECISION %s\n' % precision)  # HTVS SP XP
        if calc_rmsd == True:  # 计算与原配体rmsd 默认False
            input_file.write('CALC_INPUT_RMS True\n')
        if precision == 'XP':
            input_file.write('WRITE_XP_DESC False\n')
            input_file.write('POSTDOCK_XP_DELE 0.5\n')

    launch('glide %s_glide_dock_%s.in -JOBNAME %s-Glide-Dock-%s' %
           (pdb_code, precision, pdb_code, precision))
    os.system('mv %s-Glide-Dock-%s_pv.maegz %s_glide_dock_%s.maegz' %
              (pdb_code, precision, pdb_code, precision))
    print('\nDocking Result File:', '%s_glide_dock_%s.maegz Saved.\n' %
          (pdb_code, precision))


def dock_one_to_n(pdb_code, ligand_file, precision):
    '''
    一对多 外源配体自动对接 for multidock

    Parameters
    ----------
    pdb_code: 要对接到的受体PDB ID字符串
    ligand_file: 要对接的外源配体文件PATH
    precision: 对接精度

    '''

    ligname = ligand_file.strip().split('.')[0]
    convert_pdb_to_mae(ligand_file)
    os.chdir(pdb_code)
    grid_file = '%s_glide_grid.zip' % pdb_code
    lig_file = '../' + ligname + '.mae'

    print('\nPDB ID:', pdb_code, end='\n')
    print('Ligand Name:', ligname)
    print('Grid File:', grid_file)

    print('Prepare to Docking...\n')

    with open('%s_glide_dock_%s.in' % (ligname, precision), 'w') as input_file:
        input_file.write('GRIDFILE %s\n' % grid_file)
        input_file.write('LIGANDFILE %s\n' % lig_file)
        input_file.write('PRECISION %s\n' % precision)  # HTVS SP XP
        if precision == 'XP':
            input_file.write('WRITE_XP_DESC False\n')
            input_file.write('POSTDOCK_XP_DELE 0.5\n')

    launch('glide %s_glide_dock_%s.in -JOBNAME %s-Glide-Dock-On-%s-%s' %
           (ligname, precision, ligname, pdb_code, precision))
    os.system('mv %s-Glide-Dock-On-%s-%s_pv.maegz %s_glide_dock_on_%s_%s.maegz' %
              (ligname, pdb_code, precision, ligname, pdb_code, precision))
    print('\nDocking Result File:', '%s_glide_dock_on_%s_%s.maegz Saved.\n' %
          (ligname, pdb_code, precision))

    print('%s Docking on %s Job Complete.\n' % (ligand_file, pdb_code))
    print(''.center(80, '-'), end='\n')


def extra_data(path, precision, ligand):
    '''
    从对接完成的Maestro文件中提取数据

    Parameters
    ----------
    path: 对接完成的文件PATH
    precision: 已完成的对接精度
    ligand: 参与对接的配体名称

    Return
    ----------
    字典 Properties : Values

    '''
    file = struc.StructureReader(path)
    pdb = path.split('/')[1]
    pro_st = next(file)
    lig_st = next(file)
    prop_dic = {}

    # 需要提取的Property
    prop_dic['PDB'] = pdb
    prop_dic['Ligand'] = ligand
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

    if precision == 'XP':
        prop_dic['XP_Hbond'] = lig_st.property['r_glide_XP_HBond']

    return prop_dic


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


def preprocess(pdb):
    '''
    一些预处理操作：
    下载PDB文件 | 检查PDB晶体类型 Apo/单体/多聚体 | 是否保留单链

    Parameters
    ----------
    pdb: PDB ID字符串

    Returen
    ----------
    选择保留单链得到的文件PATH
    '''

    pdb_file = pdb + '.pdb'
    if not os.path.exists(pdb_file):  # 下载PDB文件
        getpdb.get_pdb(pdb)

    lig_lis = os.popen(
        "cat %s.pdb | grep -w -E ^HET | awk '{print $2}'" % pdb).readlines()

    if len(lig_lis) == 0:
        print('%s为Apo蛋白晶体 无法自动处理.' % pdb)
        sys.exit(1)
    elif len(lig_lis) > 1:
        print('\n')
        print('Ligand   Name  Chain Residue')
        os.system('cat %s.pdb | grep -w -E ^HET' % pdb)
        print('存在多个配体小分子 是否需要保留单链？(Y/N)')
        k = input().strip().upper()

        if k == 'Y':
            file = pdb + '.pdb'
            chain = input('输入保留链名:').strip().upper()
            pdb_file = keep_chain(pdb, file, chain)
            return pdb_file
        elif k == 'N':
            print('\nChoosed Intact Crystal.\n')
            return pdb_file
    else:
        return pdb_file


def main():
    '''
    带可控制UI的主函数
    '''

    pdb = get_pdbid()
    pdb_file = preprocess(pdb)
    lig_name = get_ligname(pdb)
    minimized_file = '%s_minimized.mae' % pdb
    grid_file = '%s_glide_grid.zip' % pdb

    print_title()
    print('\nPDB ID:', pdb, end='\n')
    print('Entry Ligand:', lig_name, end='\n')

    flag = get_flag()

    if flag == '0':
        sys.exit(0)

    elif flag == '1':
        minimized(pdb, pdb_file)
        split_com(pdb, lig_name, minimized_file)

    elif flag == '2':

        minimized_file = minimized(pdb, pdb_file)
        split_com(pdb, lig_name, minimized_file)
        grid_generate(pdb, lig_name, minimized_file)

    elif flag == '3':

        minimized_file = minimized(pdb, pdb_file)
        size = int(input('请输入Grid Box大小(Default 20Å):'))
        grid_generate(pdb, lig_name, minimized_file, size)

    elif flag == '4':

        minimized_file = minimized(pdb, pdb_file)
        lig_file = split_com(pdb, lig_name, minimized_file)[0]
        print('Ligand File:', lig_file)
        if not os.path.exists(grid_file):
            grid_file = grid_generate(pdb, lig_name, minimized_file)
        print('Grid File:', grid_file)

        dock(pdb, lig_file, grid_file, 'SP', True)

    elif flag == '5':

        minimized_file = minimized(pdb, pdb_file)
        lig_file = split_com(pdb, lig_name, minimized_file)[0]
        print('Ligand File:', lig_file)
        if not os.path.exists(grid_file):
            grid_file = grid_generate(pdb, lig_name, minimized_file)
        print('Grid File:', grid_file)
        dock(pdb, lig_file, grid_file, 'XP', True)

    elif flag == '6':

        lig_file = input('请输入配体文件PATH(单一文件 可包含多个小分子):')
        precision = input('请输入对接精度(HTVS|SP|XP):')
        convert_pdb_to_mae(lig_file)

        if not os.path.exists(grid_file):
            grid_file = grid_generate(pdb, lig_name, minimized_file)
        dock(pdb, lig_file.split('.')[0] + '.mae', grid_file, precision)

    print(''.center(80, '-'), end='\n')


def autodock(pdb, lig_name, precision):
    '''
    自动化内源配体对接 for multidock

    Parameters
    ----------
    pdb: PDB ID字符串
    lig_name: 已核对的唯一配体
    precision: 对接精度

    '''
    os.chdir(pdb)
    os.system('cp ../%s.pdb ./' % pdb)
    minimized_file = '%s_minimized.mae' % pdb
    grid_file = '%s_glide_grid.zip' % pdb

    cmd = '''cat %s.pdb | grep -w -E ^HET | awk '{if($2==\"%s\"){print $3}}' ''' % (
        pdb, lig_name)
    cmd_run = os.popen(cmd).readlines()
    if cmd_run:
        chain = re.match('[A-Z]',cmd_run[0].strip()).group()
    else:
        raise ValueError('No Match Ligand for %s' % lig_name)

    if not chain:
        raise ValueError('No Chain Match')

    pdb_file = keep_chain(pdb, pdb + '.pdb', chain)

    print('\nPDB ID:', pdb, end='\n')
    print('Entry Ligand:', lig_name, end='\n')
    print('Chain: %s' % chain)

    if not os.path.exists(minimized_file):
        minimized_file = minimized(pdb, pdb_file)

    lig_file = split_com(pdb, lig_name, minimized_file)[0]
    print('\nLigand File:', lig_file)
    if not os.path.exists(grid_file):
        grid_file = grid_generate(pdb, lig_name, minimized_file)
    print('Grid File:', grid_file)
    dock(pdb, lig_file, grid_file, precision, True)
    print('%s Self-Docking Job Complete.\n' % pdb)
    print(''.center(80, '-'), end='\n')
    return 1


def multidock(argv):
    '''
    自动多进程处理多个PDB晶体并完成自动对接 提取对接结果数据并保存为CSV文件

    Parameters
    ----------
    pdb_list: 包含多个pdb晶体ID的列表

    '''

    ligand_file = ''

    try:
        opts, argvs = getopt.getopt(argv, '-hr:l:')
    except getopt.GetoptError:
        print("usage: run py4schrodinger -r <receptors list file> -l <ligand file>\nUse run py4schrodinger.py -h for more information.")
        sys.exit(2)

    for opt, arg in opts:

        if opt == '-h':
            print(
                'usage: run py4schrodinger -r <receptors list file> [-l <ligand file>]')
            print('''
    Options
    ----------
    -r 包含所有需要对接的受体PDB ID的列表文件(.txt)
    -l 需要对接的外源配体文件(.pdb)
    
    *当没有-l参数传入时，脚本仅执行受体列表内的晶体获取、处理与内源配体对接，而不会将任何外源配体对接到列表中的受体

            ''')
            sys.exit(1)

        elif opt == '-r':
            list_file = arg
        elif opt == '-l':
            ligand_file = arg

    try:
        list_filename = list_file.split('_')[0].split('/')[-1]
        with open(list_file, 'r') as f:
            pdbs_withlig = f.readlines()
    except FileNotFoundError:
        print('Error: 未找到列表文件!')
        raise

    pdb_list = []
    for i in pdbs_withlig:
        pdb = i.split(',')[0].strip().upper()
        lig = i.split(',')[1].strip().upper()
        pdb_list.append((pdb, lig))

    print_title()
    print('\nProcessing Input List...\n')
    print('All Crystals to be Processed:',
          ''.join(str(x)+' ' for x in pdb_list))
    flag = input('\nContinue to Dock ? (Y/N)\n').strip().upper()

    if flag == 'Y':
        precision = input('请输入对接精度(HTVS|SP|XP):').strip().upper()
        print('\nNumber of Total CPU:', multiprocessing.cpu_count(), end='\n')
        cpus = int(input('请输入需要使用的CPU核心数量:'))
        multiprocessing.set_start_method('spawn')
        pool1 = multiprocessing.Pool(cpus,maxtasksperchild=1)   #每个进程必须仅使用单线程
        pool2 = multiprocessing.Pool(cpus,maxtasksperchild=1)
        dic = {}

        for pdb, lig in pdb_list:   #对接前检查

            try:
                os.makedirs(pdb)
            except FileExistsError:
                pass

            if lig == 'APO':
                dic[pdb] = lig
                continue
            if not os.path.exists(pdb + '.pdb'):  # 下载PDB文件
                getpdb.get_pdb(pdb)

            if not lig:
                dic[pdb] = get_ligname(pdb)
            else:
                dic[pdb] = lig

        total = len(dic)
        print('TOTAL:', total)

        global now_complete
        now_complete = 0

        def success_handler(result):
            global now_complete
            now_complete += result
            print("%s Job(s) Successd / Total: %s" % (now_complete,total))

        for pdb_code, ligand in dic.items():  # 采用进程池控制多线程运行
            pool1.apply_async(autodock, (pdb_code, ligand, precision,),callback=success_handler,error_callback=error_handler)

        pool1.close()  # 进程池关闭 不再提交新任务
        pool1.join()  # 阻塞进程 等待全部子进程结束

        items = dic.items()
        notpass = []
        for k,v in items:    #异常晶体跳过: APO & 共价键结合晶体
            if not os.path.exists('./%s/%s_glide_dock_%s.maegz' % (k,k,precision)):
                notpass.append((k,v))
        if notpass:
            for not_exist,lig in notpass:
                del dic[not_exist]

        if ligand_file:
            for pdb_code in dic.keys():
                pool2.apply_async(
                    dock_one_to_n, (pdb_code, ligand_file, precision,))
                #time.sleep(1.5)

            pool2.close()
            pool2.join()

        prop = ['PDB', 'Ligand', 'Docking_Score', 'rmsd', 'precision', 'ligand_efficiency', 'XP_Hbond', 'rotatable_bonds',
                'ecoul','evdw', 'emodel', 'energy', 'einternal', 'lipo', 'hbond', 'metal', 'rewards', 'erotb', 'esite']
        data = []

        for pdb, lig in dic.items():

            pro_ligand = lig
            prop_dic = extra_data('./%s/%s_glide_dock_%s.maegz' %
                                  (pdb, pdb, precision), precision, pro_ligand)
            data.append(prop_dic)
            if ligand_file:
                ligname = ligand_file.strip().split('.')[0]
                ex_dic = extra_data('./%s/%s_glide_dock_on_%s_%s.maegz' %
                                    (pdb, ligname, pdb, precision), precision, ligname)
                data.append(ex_dic)

        with open(list_filename + '_FINAL_RESULTS.csv', 'w', encoding='UTF-8', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=prop)
            writer.writeheader()
            writer.writerows(data)

        if notpass:
            print('Abandoned Crystal(s): ', notpass)
            with open('./pdbid/abandon.txt','a') as f:
                f.write('\n' + list_filename + '\n')
                for p,l in notpass:
                    f.write(p + ',' + l + '\n')
        print('\nAll Docking Jobs Done.\n')

    else:
        sys.exit(1)


if __name__ == '__main__':

    if len(sys.argv) == 1:  # UI模式
        main()

    else:  # 自动处理多晶体模式
        multidock(sys.argv[1:])
