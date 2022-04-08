import os
import logging

from pyCADD.Dock.common import launch, PDBFile, MaestroFile, GridFile

logger = logging.getLogger(__name__)

def keep_chain(pdbfile:PDBFile, chain_name:str) -> PDBFile:
    '''
    读取PDB晶体文件并将单一链的结构输出为pdb文件

    Parameters
    ----------
    pdbfile : PDBFile
        PDB晶体文件
    chain_name : str
        要保留的链名称

    Return
    ----------
    PDBFile
        保留单链结构的PDB文件
    '''

    logger.debug('Keep the single chain structure: %s' % singlechain_file)
    st = MaestroFile.get_first_structure(pdbfile.file_path)  # 读取原始PDB结构
    st_chain_only = st.chain[chain_name].extractStructure()
    singlechain_file = '%s_chain_%s.mae' % (pdbfile.pdbid, chain_name)
    st_chain_only.write(singlechain_file)
    return PDBFile(singlechain_file)

def minimize(pdbfile:PDBFile, side_chain:bool=True, missing_loop:bool=True, del_water:bool=True, overwrite:bool=False) -> MaestroFile:
    '''
    调用prepwizard模块优化PDB结构

    Parameters
    ----------
    pdbfile : PDBFile
        需要Minimize的PDB文件
    side_chain : bool
        是否优化侧链
    missing_loop : bool
        是否优化缺失的loop
    del_water : bool
        是否删除水分子
    overwrite : bool
        是否覆盖已有文件
    
    Return
    ----------
    MaestroFile
        完成优化后的Maestro文件

    '''
    
    logger.debug('Prepare to minimize %s' % pdbfile.file_name)
    pdbid = pdbfile.pdbid if pdbfile.pdbid else 'unknown'
    minimized_file = pdbid + '_minimized.mae'

    if not overwrite and os.path.exists(minimized_file):  # 如果已经进行过优化 为提高效率而跳过优化步骤
        logger.debug('File %s is existed.\n' % minimized_file)
        return MaestroFile(minimized_file)

    _job_name = '%s-Minimize' % pdbfile.pdbid
    prepwizard_command = 'prepwizard -f 3 -r 0.3 -propka_pH 7.0 -disulfides -s -j %s' % _job_name
    if side_chain:
        prepwizard_command += ' -fillsidechains'
    if missing_loop:
        prepwizard_command += ' -fillloops'
    if del_water:
        prepwizard_command += ' -watdist 0.0'

    prepwizard_command += ' %s' % pdbfile.file_path # 将pdb文件传入prepwizard
    prepwizard_command += ' %s' % minimized_file    # 将优化后的文件保存到minimized_file
    
    launch(prepwizard_command)   # 执行Minimize Job 阻塞至任务结束

    # 判断Minimized任务是否完成(是否生成Minimized结束的结构文件)
    # 无法被优化的晶体结构
    if not os.path.exists(minimized_file):  
        raise RuntimeError('%s Crystal Minimization Process Failed.' % pdbfile.pdbid)
    else:
        logger.debug('PDB minimized file: %s Saved.\n' % minimized_file)
        return MaestroFile(minimized_file)

def grid_generate(maestrofile:MaestroFile, ligname:str, gridbox_size:int=20, overwrite:bool=False) -> GridFile:
    '''
    自动编写glide grid输入文件并启动Glide Grid生成任务

    Parameters
    ----------
    maestrofile : MaestroFile
        需要生成Grid的Maestro文件
    lig_name : str
        作为坐标参考的Ligand名称
    gridbox_size : int
        grid box大小 默认20Å
    overwrite : bool
        是否覆盖已有文件

    Return
    ----------
    GridFile
        生成的格点文件

    '''

    lig_molnum = maestrofile.get_lig_molnum(ligname)
    pdbid = maestrofile.pdbid if maestrofile.pdbid else 'unknown'
    grid_file = f'{pdbid}_glide-grid_{ligname}.zip'
    logger.debug(f'Prepare to generate grid file: {grid_file}')

    # 如果已经生成了格点文件则跳过生成过程
    if os.path.exists(grid_file) and not overwrite: 
        logger.debug('File %s is existed.\n' % grid_file)
        return GridFile(grid_file)

    input_file = f'{pdbid}_glide-grid_{ligname}.in'
    job_name = f'{pdbid}-{ligname}-Grid-Generate'
    outsize = gridbox_size + 10

    # 编写glide输入文件
    glide_grid_config = [
        'GRIDFILE %s\n' % grid_file,
        'INNERBOX 10,10,10\n',
        'OUTERBOX %d,%d,%d \n' % (outsize, outsize, outsize),
        'LIGAND_MOLECULE %s\n' % lig_molnum,
        'RECEP_FILE %s\n' % maestrofile.file_path
    ]

    with open(input_file, 'w') as f:  
        f.writelines(glide_grid_config)

    launch(f'glide {input_file} -JOBNAME {job_name}')

    logger.debug('Grid File %s Generated.\n' % grid_file)
    return GridFile(grid_file)

def dock(lig_file:MaestroFile, grid_file:GridFile, precision:str='SP', calc_rmsd:bool=False, overwrite:bool=False) -> MaestroFile:
    '''
    一对一/多对一 glide dock任务输入文件编写与运行

    Parameters
    ----------
    lig_file : MaestroFile
        配体文件
    grid_file_path : GridFile
        格点文件
    precision : str
        对接精度(HTVS|SP|XP) 默认SP
    calc_rmsd : bool
        是否计算rmsd to input ligand geometries 默认False
    overwrite : bool
        是否覆盖已有文件

    Return
    ---------
    MaestroFile
        对接结果文件
    '''
    pdbid = grid_file.pdbid if grid_file.pdbid else 'unknown'
    # 共结晶配体名称
    internal_ligand = grid_file.internal_ligand
    # 正在对接的配体名称
    docking_ligand = lig_file.file_prefix
    dock_result_file = f'{pdbid}_{internal_ligand}_glide-dock_{docking_ligand}_{precision}.maegz'
    
    logger.debug(f'Prepare to dock {docking_ligand} on {pdbid}')
    # 如果已有对接成功文件 跳过对接步骤
    if os.path.exists(dock_result_file) and not overwrite:     
        logger.debug('File %s is existed.\n' % dock_result_file)                           
        return MaestroFile(dock_result_file)

    input_file = f'{pdbid}_{internal_ligand}_glide-dock_{docking_ligand}_{precision}.in'
    job_name = f'{pdbid}-{internal_ligand}-Glide-Dock-{docking_ligand}-{precision}'
    output_file = job_name + '_pv.maegz'
    glide_dock_config = [
        'GRIDFILE %s\n' % grid_file.file_path,
        'LIGANDFILE %s\n' % lig_file.file_path,
        'PRECISION %s\n' % precision
    ]
    if calc_rmsd is True:
        glide_dock_config.append('CALC_INPUT_RMS True\n')
    if precision == 'XP':
        glide_dock_config.append('WRITE_XP_DESC False\n')
        glide_dock_config.append('POSTDOCK_XP_DELE 0.5\n')

    with open(input_file, 'w') as f:
        f.writelines(glide_dock_config)

    launch(f'glide {input_file} -JOBNAME {job_name}')

    try:
        os.rename(output_file, dock_result_file)
    except Exception:
        logger.debug(f'{pdbid}-{docking_ligand} Docking Failed')
        return None

    logger.debug(f'{pdbid}-{docking_ligand} Glide Docking Completed')
    logger.debug(f'Docking Result File: {dock_result_file} Saved.')

    return MaestroFile(dock_result_file)

def calc_mmgbsa(complex_file:MaestroFile, overwrite:bool=False) -> MaestroFile:
    '''
    计算MM-GBSA结合能

    Parameters
    ----------
    complex_file : MaestroFile
        需要计算MMGBSA的复合物文件

    Return
    ----------
    Maestro
        计算MM-GB/SA完成的复合物文件

    '''
    
    # 按照dock()自动生成的文件名获取相关信息
    prefix = complex_file.file_prefix

    logger.debug('Prepare to Calculate MM-GB/SA Binding Energy : %s' % prefix)
    mmgbsa_result_file = prefix + '_mmgbsa.maegz'

    # 已计算则跳过计算过程
    if os.path.exists(mmgbsa_result_file) and not overwrite:
        logger.debug('File %s is existed.\n' % mmgbsa_result_file)
        return MaestroFile(mmgbsa_result_file)

    launch(f'prime_mmgbsa -j {prefix} {complex_file.file_path}')
    
    try:
        os.rename(f'{prefix}-out.magez', mmgbsa_result_file)
    except Exception:
        logger.debug(f'{prefix} Prime MM-GB/SA Calculating Failed')
        return None

    logger.debug(f'MM-GB/SA Calculating Result File: {mmgbsa_result_file} Saved.\n')

    return MaestroFile(mmgbsa_result_file)

def calc_volume(recep_file:MaestroFile, lig_file:MaestroFile, overwrite:bool=False) -> MaestroFile:
    '''
    Sitemap计算结合口袋体积

    Parameters
    ----------
    recep_file : MaestroFile
        进行口袋体积分析的受体文件
    lig_file : MaestroFile
        定义口袋位置的配体文件
        
    Return
    ----------
    MaestroFile
        sitemap计算完成的文件
    '''
    pdbid = recep_file.pdbid if recep_file.pdbid else 'unknown'
    sitemap_result_file = '%s_sitemap_out.maegz' % pdbid
    logger.debug('Prepare to calculate volume of site in %s' % pdbid)

    if os.path.exists(sitemap_result_file) and not overwrite:
        logger.debug('File %s is existed.\n' % sitemap_result_file)
        return MaestroFile(sitemap_result_file)

    launch(f'sitemap -sitebox 6 -keeplogs yes -ligmae {lig_file.file_path} -prot {recep_file.file_path} -j {pdbid}_sitemap')
    
    if not os.path.exists(sitemap_result_file):
        logger.debug(f'{pdbid} Sitemap Calculating Failed')
        return None

    logger.debug('Sitemap Calculating File: %s Saved.\n' % sitemap_result_file)
    return MaestroFile(sitemap_result_file)

def calc_admet(lig_file:MaestroFile, overwrite:bool=False) -> MaestroFile:
    '''
    计算单个化合物/配体的ADMET特征描述符

    Parameter  
    ----------
    lig_file : MaestroFile
        化合物/配体文件

    Return
    ----------
    MaestroFile
        ADMET 计算结果文件
    '''
    prefix = lig_file.file_prefix + '_ADMET'
    admet_result_file = prefix + '.mae'

    logger.debug('Prepare to calculate ADMET properties: %s' % prefix)

    if os.path.exists(admet_result_file) and not overwrite:
        logger.debug('File %s is existed.\n' % admet_result_file)
        return MaestroFile(admet_result_file)

    launch(f'qikprop -outname {prefix} {lig_file.file_path}')
    try:
        os.rename(prefix + '.mae', admet_result_file)
    except Exception:
        logger.debug(f'{prefix} ADMET Calculating Failed')
        return None

    logger.debug('Qikprop Calculation File: %s Saved.\n' % admet_result_file)
    return MaestroFile(admet_result_file)
    