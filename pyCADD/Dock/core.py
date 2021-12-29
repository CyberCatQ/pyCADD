import os
import logging

from pyCADD.Dock.prepare import convert_format, load_st
from pyCADD.utils.check import check_file, checkpdb
from pyCADD.utils.getinfo import get_ligmol_info
from schrodinger.application.glide import poseviewconvert as pvc
from schrodinger.job import jobcontrol as jc

logger = logging.getLogger('pyCADD.Dock.core')

def launch(cmd:str, timeout:int=None):
    '''
        使用jobcontrol启动一项job并等待结束

        Parameters
        ----------
        cmd : str
            等待执行的命令字符串
        timeout : int
            超时时限(秒)
        '''

    cmd_list = cmd.split(' ')  # launch_job以列表形式提交参数
    job = jc.launch_job(cmd_list, timeout=timeout)
    logger.debug('Command: %s' % cmd)
    logger.debug('JobId: %s' % job.JobId)
    logger.debug('Job Name: %s' % job.Name)
    logger.debug('Job Status: %s\n' % job.Status)
    job.wait()  # 阻塞进程 等待Job结束

    # 如果任务失败
    if not job.StructureOutputFile:
        logger.debug('Job %s Failed' % job.Name)
        return
    else:
        logger.debug('File %s saved.' % job.StructureOutputFile)


def keep_chain(pdbfile_path:str, chain_name:str) -> str:
    '''
    读取PDB晶体文件并将单一链的结构输出为pdb文件

    Parameters
    ----------
    pdbfile : str
        待处理原始文件PATH
    chain_name : str
        要保留的链名称

    Return
    ----------
    str
        保留单链结构的文件名
    '''
    pdbfile = os.path.basename(pdbfile_path)
    st = load_st(pdbfile_path)  # 读取原始PDB结构
    st_chain_only = st.chain[chain_name].extractStructure()
    singlechain_file = '%s_chain_%s.mae' % (pdbfile.split('.')[0], chain_name)
    st_chain_only.write(singlechain_file)
    logger.debug('Keep the single chain structure: %s' % singlechain_file)
    return singlechain_file

def minimize(pdbfile_path:str) -> str:
    '''
    调用prepwizard模块自动优化PDB结构文件

    Parameters
    ----------
    pdbfile_path : str
        需要优化的文件PATH
    
    Return
    ----------
    str
        完成优化后的文件名

    '''
    
    pdbfile = os.path.basename(pdbfile_path)
    logger.debug('Prepare to minimize %s' % pdbfile)
    minimized_file = str(pdbfile.split('.')[0] + '_minimized.mae')

    if os.path.exists(minimized_file):  # 如果已经进行过优化 为提高效率而跳过优化步骤
        logger.debug('File %s is existed.\n' % minimized_file)
        return minimized_file

    prepwizard_command = 'prepwizard -f 3 -r 0.3 -propka_pH 7.0 -disulfides -s -j %s-Minimize %s %s' % (pdbfile.split('.')[0],
                                                                                                            pdbfile_path, minimized_file)
    launch(prepwizard_command)   # 阻塞至任务结束

    # 判断Minimized任务是否完成(是否生成Minimized结束的结构文件)
    # 无法被优化的晶体结构
    if not os.path.exists(minimized_file):  
        raise RuntimeError('%s Crystal Minimization Process Failed.' % pdbfile.split('.')[0])  
    else:
        logger.debug('PDB minimized file: %s Saved.\n' % minimized_file)
        return minimized_file

def grid_generate(pdbid:str, ligname:str, st_file_path:str, gridbox_size:int=20) -> str:
    '''
    自动编写glide grid输入文件并启动Glide Grid生成任务

    Parameters
    ----------
    pdbid : str
        PDB ID
    ligname : str
        配体ID(RCSB ID)
    st_file_path : str
        需要构建grid box的文件PATH
    gridbox_size : int
        grid box大小 默认20Å

    Return
    ----------
    str
        生成的格点文件名

    '''
        
    lig_molnum = get_ligmol_info(st_file_path, ligname)
    outsize = gridbox_size + 10

    grid_file = pdbid + '_glide_grid_%s.zip' % ligname
    logger.debug('Prepare to generate grid file: %s' % grid_file)

    # 如果已经生成了格点文件则跳过生成过程
    if os.path.exists(grid_file):   
        logger.debug('File %s is existed.\n' % grid_file)
        return grid_file

    # 编写glide输入文件
    with open('%s_grid_generate_%s.in' % (pdbid, ligname), 'w') as input_file:  
        input_file.write('GRIDFILE %s\n' % grid_file)           # 输出的Grid文件名
        input_file.write('INNERBOX 10,10,10\n')                 # Box大小参数
        input_file.write('OUTERBOX %d,%d,%d \n' %
                            (outsize, outsize, outsize))        # Box大小参数
        input_file.write('LIGAND_MOLECULE %s\n' %
                            lig_molnum)                         # 识别Ligand并设定grid box中心为质心
        input_file.write('RECEP_FILE %s' % st_file_path)        # structure输入文件PATH
    launch('glide %s_grid_generate_%s.in -JOBNAME %s-%s-Grid-Generate' % (pdbid, ligname, pdbid, ligname))

    logger.debug('Grid File %s Generated.\n' % grid_file)
    return grid_file

def split_com(complex_file_path:str, ligname:str) -> tuple:
    '''
    拆分复合体为配体和受体

    Parameters
    ----------
    complex_file_path : str
        待拆分的复合物文件PATH
    ligname : str
        配体ID(RCSB ID)

    Return
    ----------
    tuple
        (ligfile_PATH, recepfile_PATH)

        '''

    complex_file = os.path.basename(complex_file_path)
    pdbid = complex_file.split('_')[0]

    if not checkpdb(pdbid):
        raise ValueError('Invaild PDB ID : %s' % pdbid)

    st = load_st(complex_file_path)

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
        logger.debug('There are %s "%s" in %s: %s' % (len(residue_list), ligname, complex_file, resname_list))
        logger.debug('The First One is Selected.\n')

    lig_file = '%s-lig-%s.mae' % (pdbid, ligname)
    recep_file = '%s-pro-%s.mae' % (pdbid, ligname)     
    com_file = '%s-com-%s.pdb' % (pdbid, ligname)
    # 生成并保存配体独立mae文件
    comp.writeLigand(lig_file) 
    lig_st = load_st(lig_file)    
    lig_st.property['b_user_Activity'] = 1
    lig_st.write(lig_file)
    logger.debug('Ligand file %s is generated.' % lig_file)     
    comp.writeReceptor(recep_file)
    logger.debug('Receptor file %s is generated.' % recep_file)
    st.write(com_file)
    logger.debug('Complex file %s is generated.' % com_file)

    # 自动生成PDB格式
    convert_format(lig_file, 'pdb')
    convert_format(recep_file, 'pdb')
    logger.debug('PDB format converted.')

    return lig_file, recep_file

def dock(lig_file_path:str, grid_file_path:str, precision:str='SP', calc_rmsd:bool=False) -> str:
    '''
    一对一/多对一 glide dock任务输入文件编写与运行

    Parameters
    ----------
    lig_file_path : str
        配体文件PATH
    grid_file_path : str
        格点文件PATH
    precision : str
        对接精度(HTVS|SP|XP) 默认SP
    calc_rmsd : bool
        是否计算rmsd to input ligand geometries 默认False

    Return
    ---------
    str
        对接结果文件名
    '''
    grid_file_name = os.path.basename(grid_file_path).split('.')[0]
    pdbid = grid_file_name.split('_')[0]
    # 共结晶配体名
    internal_ligand = grid_file_name.split('_')[3]

    # 正在对接的配体名称
    lig_name = os.path.basename(lig_file_path).split('.')[0].split('_')[-1]
    logger.debug('Prepare to dock %s on %s' % (lig_name, pdbid))

    dock_file = '%s_%s_glide_dock_%s_%s.maegz' % (pdbid, internal_ligand, lig_name, precision)
    # 如果已有对接成功文件 跳过对接步骤
    if os.path.exists(dock_file):     
        logger.debug('File %s is existed.\n' % dock_file)    
        logger.info('%s-%s Glide Docking Completed' % (pdbid, lig_name))                         
        return dock_file

    with open('%s_%s_glide_dock_%s_%s.in' % (pdbid, internal_ligand, lig_name, precision), 'w') as input_file:
        input_file.write('GRIDFILE %s\n' % grid_file_path)
        # 如果lig文件包含多个lig mol 将自动从一至末尾全部dock
        input_file.write('LIGANDFILE %s\n' % lig_file_path)
        # 对接精度 HTVS SP XP
        input_file.write('PRECISION %s\n' % precision)  
        # 计算与原配体rmsd 默认False
        if calc_rmsd == True:                           
            input_file.write('CALC_INPUT_RMS True\n')
        if precision == 'XP':
            input_file.write('WRITE_XP_DESC False\n')
            input_file.write('POSTDOCK_XP_DELE 0.5\n')

    launch('glide %s_%s_glide_dock_%s_%s.in -JOBNAME %s-%s-Glide-Dock-%s-%s' % (pdbid, internal_ligand, lig_name, precision, pdbid, internal_ligand, lig_name, precision))
    if not check_file('%s-%s-Glide-Dock-%s-%s_pv.maegz' % (pdbid, internal_ligand, lig_name, precision)):
        logger.warning('%s-%s Glide Docking Failed' % (pdbid, lig_name))
        return
    os.system('mv %s-%s-Glide-Dock-%s-%s_pv.maegz %s_%s_glide_dock_%s_%s.maegz' % (pdbid, internal_ligand, lig_name, precision, pdbid, internal_ligand, lig_name, precision))
    logger.info('%s-%s Glide Docking Completed' % (pdbid, lig_name))
    logger.debug('Docking Result File: %s_%s_glide_dock_%s_%s.maegz Saved.' % (pdbid, internal_ligand, lig_name, precision))

    return dock_file

def cal_mmgbsa(dock_file_path:str) -> str:
    '''
    计算MM-GBSA结合能

    Parameters
    ----------
    dock_file_path : str
        对接完成的结果文件PATH(包含受体与配体)   

    Return
    ----------
    str
        计算MM-GB/SA完成的复合物文件名

    '''
    dock_file_name = os.path.basename(dock_file_path).split('.')[0]

    # 按照dock()自动生成的文件名获取相关信息
    pdbid = dock_file_name.split('_')[0]
    ligname = dock_file_name.split('_')[4]

    logger.debug('Prepare to Calculate MM-GB/SA Binding Energy : %s-%s' % (pdbid, ligname))
    mmgbsa_file = dock_file_name + '_mmgbsa.maegz'

    # 已计算则跳过计算过程
    if os.path.exists(mmgbsa_file):
        logger.debug('File %s is existed.\n' % mmgbsa_file)
        return mmgbsa_file

    launch('prime_mmgbsa -j %s %s' % (mmgbsa_file.split('.')[0], dock_file_path))
    cmd = os.system('mv %s-out.maegz %s'% (mmgbsa_file.split('.')[0], mmgbsa_file))
    if cmd != 0:
        raise RuntimeError('%s-%s Prime MM-GB/SA Calculating Failed.' % (pdbid, ligname))
        
    logger.debug('MM-GB/SA Calculating Result File: %s Saved.\n' % mmgbsa_file)

    return mmgbsa_file

def cal_volume(recep_file_path:str, lig_file_path:str) -> str:
    '''
    Sitemap计算结合口袋体积

    Parameters
    ----------
    recep_file_path : str
        进行口袋体积分析的受体文件PATH
    lig_file_path : str
        定义口袋位置的配体文件PATH
        
    Return
    ----------
    str
        sitemap计算完成的文件名
    '''
    pdbid = os.path.basename(recep_file_path).split('-')[0]
    sitemap_file = '%s_sitemap_out.maegz' % pdbid
    logger.debug('Prepare to calculate volume of site in %s' % pdbid)
    if os.path.exists(sitemap_file):
        logger.debug('File %s is existed.\n' % sitemap_file)
        return sitemap_file

    launch('sitemap -sitebox 6 -keeplogs yes -ligmae %s -prot %s -j %s_sitemap' % (lig_file_path, recep_file_path, pdbid))
    
    if not os.path.exists(sitemap_file):
        raise RuntimeError('%s Sitemap Calculating Failed.' % pdbid)

    logger.debug('Sitemap Calculating File: %s Saved.\n' % sitemap_file)
    return sitemap_file

def cal_admet(lig_file_path:str):
    '''
    计算单个化合物/配体的ADMET特征描述符

    Parameter  
    ----------
    lig_file_path : str
        化合物/配体文件PATH

    Return
    ----------
    str
        ADMET 计算结果文件PATH
    '''

    lig_file = os.path.basename(lig_file_path)
    prefix = lig_file.split('.')[0] + '_ADMET'

    admet_file = prefix + '.mae'
    logger.debug('Prepare to calculate ADMET properties: %s' % lig_file)
    if os.path.exists(admet_file):
        logger.debug('File %s is existed.\n' % admet_file)
        return admet_file

    logger.debug('Calculating ADMET Propeties of %s' % lig_file)

    launch('qikprop -outname %s %s' % (prefix, lig_file_path))
    os.system('mv %s-out.mae %s' % (prefix, admet_file))

    logger.debug('Qikprop Calculation File: %s Saved.\n' % admet_file)

    return admet_file
    