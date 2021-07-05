import object_py4schrodinger


class mutilconsole:

    def __dock_one_to_n(self, origin_ligname, ex_ligand_file, precision):
        '''
        一对多 外源配体自动对接 for multidock

        Parameters
        ----------
        pdb_code: 要对接到的受体PDB ID字符串
        origin_ligname: 对接位置原配体名称
        ex_ligand_file: 要对接的外源配体文件PATH
        precision: 对接精度

        '''

        pdbid = self.pdbid
        origin_ligname = self.ligname
        ex_ligname = ex_ligand_file.strip().split('.')[0]
        self.convert_format(ex_ligand_file, 'mae')

        os.chdir('./database/' + pdbid)
        grid_file = '%s_glide_grid_%s.zip' % (pdbid, origin_ligname)
        lig_file = '../../' + ligname + '.mae'

        print('\nPDB ID:', pdb_code, end='\n')
        print('Ligand Name:', ligname)
        print('Grid File:', grid_file)

        print('Prepare to Docking...\n')

        with open('%s_glide_dock_%s_%s.in' % (ligname, origin_ligname, precision), 'w') as input_file:
            input_file.write('GRIDFILE %s\n' % grid_file)
            input_file.write('LIGANDFILE %s\n' % lig_file)
            input_file.write('PRECISION %s\n' % precision)  # HTVS SP XP
            if precision == 'XP':
                input_file.write('WRITE_XP_DESC False\n')
                input_file.write('POSTDOCK_XP_DELE 0.5\n')

        launch('glide %s_glide_dock_%s_%s.in -JOBNAME %s-Glide-Dock-On-%s-%s-%s' %
            (ligname, origin_ligname, precision, ligname, pdb_code, origin_ligname, precision))

        if not os.path.exists('%s-Glide-Dock-On-%s-%s-%s_pv.maegz' % (ligname, pdb_code, origin_ligname, precision)):
            raise RuntimeError('%s-%s-%s Gilde Dock Failed' % (ligname, pdb_code, origin_ligname))

        os.system('mv %s-Glide-Dock-On-%s-%s-%s_pv.maegz %s_glide_dock_on_%s_%s_%s.maegz' %
                (ligname, pdb_code, origin_ligname, precision, ligname, pdb_code, origin_ligname, precision))
        print('\nDocking Result File:', '%s_glide_dock_on_%s_%s_%s.maegz Saved.\n' %
            (ligname, pdb_code, origin_ligname, precision))

        print('%s Docking on %s-%s Job Complete.\n' %
            (ligname, pdb_code, origin_ligname))
        print(''.center(80, '-'), end='\n')
        return 1

def autodock(pdb, lig_name, precision):
    '''
    自动化内源配体对接 for multidock

    Parameters
    ----------
    pdb: PDB ID字符串
    lig_name: 已核对的唯一配体
    precision: 对接精度

    '''
    os.chdir('./database/' + pdb)
    if not os.path.exists(pdb + '.pdb'):
        os.system('cp ../../%s.pdb ./' % pdb)

    cmd = '''cat %s.pdb | grep -w -E ^HET | awk '{if($2==\"%s\"){print $3}}' ''' % (
        pdb, lig_name)
    cmd_run = os.popen(cmd).readlines()
    if cmd_run:
        chain = re.match('[A-Z]', cmd_run[0].strip()).group()
    else:
        raise ValueError('No Match Ligand for %s' % lig_name)

    if not chain:
        raise ValueError('No Chain Match')

    pdb_file = keep_chain(pdb, pdb + '.pdb', chain)

    print('\nPDB ID:', pdb, end='\n')
    print('Entry Ligand:', lig_name, end='\n')
    print('Chain: %s' % chain)

    minimized_file = minimized(pdb, pdb_file)
    lig_file = split_com(pdb, lig_name, minimized_file)[0]
    print('\nLigand File:', lig_file)
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
    flag = ''
    cpus = ''
    precision = ''

    try:
        opts, argvs = getopt.getopt(argv, '-hkr:l:p:n:')
    except getopt.GetoptError:
        print(
            "usage: run py4schrodinger -r <receptors list file> [-l <ligand file> -p <precision> -n <cpus> -k]\nUse run py4schrodinger.py -h for more information.")
        sys.exit(2)

    for opt, arg in opts:

        if opt == '-h':
            print(
                'usage: run py4schrodinger -r <receptors list file> [-l <ligand file> -p <precision> -n <cpus> -k]')
            print('''
    Options
    ----------
    -r 包含所有需要对接的受体PDB ID的列表文件(.txt)
    -l 需要对接的外源配体文件(.pdb)
    -p 对接精度(HTVS|SP|XP)
    -n 要使用的CPU数量
    -k 不进行对接启动确认 直接开始对接
    
    *当没有-l参数传入时，脚本仅执行受体列表内的晶体获取、处理与内源配体对接，而不会将任何外源配体对接到列表中的受体

            ''')
            sys.exit(1)

        elif opt == '-r':
            list_file = arg
        elif opt == '-l':
            ligand_file = arg
        elif opt == '-p':
            precision = arg
        elif opt == '-k':
            flag = 'Y'
        elif opt == '-n':
            cpus = int(arg)

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

    global total
    total = len(pdb_list)
    print('\nTOTAL CRYSTALS:', total)

    if not flag:
        flag = input('\nContinue to Dock ? (Y/N)\n').strip().upper()

    if flag == 'Y':
        if not precision:
            precision = input('请输入对接精度(HTVS|SP|XP):').strip().upper()

        if not cpus:
            print('\nNumber of Total CPU:',
                  multiprocessing.cpu_count(), end='\n')
            cpus = int(input('请输入需要使用的CPU核心数量:'))

        print('Using Number of CPU: %s' % cpus)
        print('Docking Precision: %s' % precision)

        multiprocessing.set_start_method('spawn')
        pool1 = multiprocessing.Pool(cpus, maxtasksperchild=1)  # 每个进程必须仅使用单线程
        pool2 = multiprocessing.Pool(cpus, maxtasksperchild=1)

        for pdb, lig in pdb_list:  # 对接前检查

            try:
                os.makedirs('./database/' + pdb)
            except FileExistsError:
                pass

            if not os.path.exists('./database/%s/%s.pdb' % (pdb, pdb)):  # 下载PDB文件
                getpdb.get_pdb(pdb)

        global now_complete
        now_complete = 0

        for pdb_code, ligand in pdb_list:  # 采用进程池控制多线程运行
            pool1.apply_async(autodock, (pdb_code, ligand, precision,),
                              callback=success_handler, error_callback=error_handler)
            time.sleep(0.5)

        pool1.close()  # 进程池关闭 不再提交新任务
        pool1.join()  # 阻塞进程 等待全部子进程结束

        items = pdb_list
        notpass = []
        for k, v in items:  # 异常晶体跳过: APO & 共价键结合晶体
            if not os.path.exists('./database/%s/%s_glide_dock_%s_%s.maegz' % (k, k, v, precision)):
                notpass.append((k, v))
        if notpass:
            for not_exist, lig in notpass:
                pdb_list.remove((not_exist, lig))

        withlig = ''
        if ligand_file:
            
            withlig = '_' + ligand_file.split('.')[0]
            print('\n')
            print(''.center(80,'-'))
            print('User-Defined Ligand Docking'.center(80))
            print(''.center(80,'-'))
            print('\nUser-Defined Ligand:', ligand_file.split('.')[0])

            now_complete = 0
            for pdb_code, lig in pdb_list:
                pool2.apply_async(dock_one_to_n, (pdb_code, lig, ligand_file, precision,),
                                  callback=success_handler, error_callback=error_handler)
                time.sleep(0.5)

            pool2.close()
            pool2.join()

        prop_xp = ['PDB', 'Ligand', 'Docking_Score', 'rmsd', 'precision', 'ligand_efficiency', 'XP_Hbond', 'rotatable_bonds',
                'ecoul', 'evdw', 'emodel', 'energy', 'einternal']
        prop_sp = ['PDB', 'Ligand', 'Docking_Score', 'rmsd', 'precision', 'ligand_efficiency', 'rotatable_bonds',
                'ecoul', 'evdw', 'emodel', 'energy', 'einternal','lipo', 'hbond', 'metal', 'rewards', 'erotb', 'esite']

        data = []
        dock_fail = []

        for pdb, lig in pdb_list:

            pro_ligand = lig
            prop_dic = extra_data('./database/%s/%s_glide_dock_%s_%s.maegz' %
                                  (pdb, pdb, pro_ligand, precision), precision, pro_ligand)
            data.append(prop_dic)
            if ligand_file:
                ligname = ligand_file.strip().split('.')[0]
                dock_result_file = './database/%s/%s_glide_dock_on_%s_%s_%s.maegz' % (
                    pdb, ligname, pdb, pro_ligand, precision)
                if os.path.exists(dock_result_file):
                    ex_dic = extra_data(dock_result_file, precision, ligname)
                    data.append(ex_dic)
                else:
                    dock_fail.append((pdb,pro_ligand,ligname))

        with open(list_filename + '_FINAL_RESULTS%s.csv' % withlig, 'w', encoding='UTF-8', newline='') as f:
            if precision == 'XP':
                writer = csv.DictWriter(f, fieldnames=prop_xp)
            elif precision == 'SP':
                writer = csv.DictWriter(f, fieldnames=prop_sp)
            writer.writeheader()
            writer.writerows(data)

        if notpass:
            print('\nAbandoned Crystal(s): ', notpass)
            with open('./pdbid/abandon.txt', 'a') as f:
                f.write('\n' + list_filename + '\n')
                for p, l in notpass:
                    f.write(p + ',' + l + '\n')
        if dock_fail:
            print('Docking Failed Crystal(s): ', dock_fail)

        print('\nAll Docking Jobs Done.\n')

    else:
        sys.exit(1)