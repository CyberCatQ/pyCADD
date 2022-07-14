import os
import logging
import shutil

from pyCADD.Dock.common import launch, PDBFile, GridFile, LigandFile, ReceptorFile, ComplexFile, DockResultFile
from pyCADD.Dock.config import GLIDE_FORCEFIELD as FORCEFIELD

logger = logging.getLogger(__name__)

def minimize(pdbfile:PDBFile, side_chain:bool=True, missing_loop:bool=True, del_water:bool=True, save_dir:str=None, overwrite:bool=False) -> ComplexFile:
    '''
    调用prepwizard模块优化PDB结构

    Parameters
    ----------
    pdbfile : PDBFile | MaestroFile
        需要Minimize的PDB|Mae文件
    side_chain : bool
        是否优化侧链
    missing_loop : bool
        是否优化缺失的loop
    del_water : bool
        是否删除水分子
    save_dir : str
        保存优化后的文件的目录
    overwrite : bool
        是否覆盖已有文件
    
    Return
    ----------
    ComplexFile
        完成优化后的Maestro文件

    '''
    

    logger.debug('Prepare to minimize %s' % pdbfile.file_name)
    # pdbid = pdbfile.pdbid if pdbfile.pdbid else 'unknown'
    ligand_id = pdbfile.ligid if pdbfile.ligid else None
    ligand_resnum = pdbfile.lig_resnum if pdbfile.lig_resnum else None
    
    save_dir = save_dir if save_dir else os.getcwd()
    minimized_file = os.path.join(save_dir, f'{pdbfile.file_prefix}_minimized.mae')

    _cwd = os.getcwd()
    os.chdir(save_dir)

    if not overwrite and os.path.exists(minimized_file):  # 如果已经进行过优化 为提高效率而跳过优化步骤
        logger.debug('File %s is existed.' % minimized_file)
        return ComplexFile(minimized_file, ligand_id, ligand_resnum)

    _job_name = '%s-Minimize' % pdbfile.file_prefix
    prepwizard_command = 'prepwizard -f 3 -r 0.3 -propka_pH 7.0 -disulfides -s -j %s' % _job_name
    if side_chain:
        prepwizard_command += ' -fillsidechains'
    if missing_loop:
        prepwizard_command += ' -fillloops'
    if del_water:
        prepwizard_command += ' -watdist 0.0'

    # Prepwizard要求输出文件名必须为相对路径 只能先输出 后移动
    _minimized_file = f'{pdbfile.file_prefix}_minimized.mae'

    prepwizard_command += ' %s' % pdbfile.file_path # 将pdb文件传入prepwizard
    prepwizard_command += ' %s' % _minimized_file    # 将优化后的文件保存到minimized_file
    
    launch(prepwizard_command)   # 执行Minimize Job 阻塞至任务结束
    
    # 判断Minimized任务是否完成(是否生成Minimized结束的结构文件)
    try:
        shutil.move(_minimized_file, minimized_file)
    except FileNotFoundError:
        raise RuntimeError('%s Crystal Minimization Process Failed.' % pdbfile.file_prefix)
    else:
        os.chdir(_cwd)
        logger.debug('PDB minimized file: %s Saved.' % minimized_file)
        return ComplexFile(minimized_file, ligand_id, ligand_resnum)

def grid_generate(complex_file:ComplexFile, gridbox_size:int=20, save_dir:str=None, overwrite:bool=False) -> GridFile:
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
    save_dir : str
        保存Grid文件的目录
    overwrite : bool
        是否覆盖已有文件

    Return
    ----------
    GridFile
        生成的格点文件

    '''
    ligname = complex_file.ligid
    lig_resnum = int(complex_file.lig_resnum) if complex_file.lig_resnum else None
    if ligname is None and lig_resnum is None:
        ligname = input('Ligand name is not specified in grid file. Please input the ligand name: ')

    lig_molnum = complex_file.get_lig_molnum(ligname, lig_resnum)
    pdbid = complex_file.pdbid if complex_file.pdbid else 'unknown'
    save_dir = save_dir if save_dir else os.getcwd()
    grid_file = os.path.join(save_dir, f'{pdbid}_glide-grid_{ligname}.zip')
    logger.debug(f'Prepare to generate grid file: {grid_file}')
    logger.debug('Grid center is set to ligand %s %s' % (ligname, lig_resnum))

    _cwd = os.getcwd()
    os.chdir(save_dir)
    # 如果已经生成了格点文件则跳过生成过程
    if os.path.exists(grid_file) and not overwrite: 
        logger.debug('File %s is existed.' % grid_file)
        return GridFile(grid_file, ligname, lig_resnum)

    input_file = f'{pdbid}_glide-grid_{ligname}.in'
    job_name = f'{pdbid}-{ligname}-Grid-Generate'
    outsize = gridbox_size + 10

    # 编写glide输入文件
    glide_grid_config = [
        'FORCEFIELD %s\n' % FORCEFIELD,
        'GRIDFILE %s\n' % grid_file,
        'INNERBOX 10,10,10\n',
        'OUTERBOX %d,%d,%d \n' % (outsize, outsize, outsize),
        'LIGAND_MOLECULE %s\n' % lig_molnum,
        'RECEP_FILE %s\n' % complex_file.file_path
    ]

    with open(input_file, 'w') as f:  
        f.writelines(glide_grid_config)

    launch(f'glide {input_file} -JOBNAME {job_name}')
    if not os.path.exists(grid_file):
        raise RuntimeError(f'{grid_file} Generation Failed.')
    else:
        logger.debug('Grid File %s Generated.' % grid_file)
        os.chdir(_cwd)

        return GridFile(grid_file, ligname, lig_resnum)

def dock(grid_file:GridFile, lig_file:LigandFile, precision:str='SP', calc_rmsd:bool=False, ligand_only:bool=False, save_dir:str=None, overwrite:bool=False) -> DockResultFile:
    '''
    一对一/多对一 glide dock任务输入文件编写与运行

    Parameters
    ----------
    grid_file_path : GridFile
        格点文件
    lig_file : LigandFile
        配体文件
    precision : str
        对接精度(HTVS|SP|XP) 默认SP
    calc_rmsd : bool
        是否计算rmsd to input ligand geometries 默认False
    save_dir : str
        保存Dock结果文件的目录 保存于该目录的PDBID文件夹下
    overwrite : bool
        是否覆盖已有文件

    Return
    ---------
    DockResultFile
        对接结果文件
    '''
    pdbid = grid_file.pdbid if grid_file.pdbid else 'unknown'
    # 共结晶配体名称
    internal_ligand = grid_file.internal_ligand
    internal_ligand_resnum = grid_file.lig_resnum
    # 正在对接的配体名称
    docking_ligand = lig_file.file_prefix
    save_dir = save_dir if save_dir else os.getcwd()
    save_dir = os.path.join(save_dir, pdbid)
    os.makedirs(save_dir, exist_ok=True)
    dock_result_file = os.path.join(save_dir, f'{pdbid}_{internal_ligand}_glide-dock_{docking_ligand}_{precision}.maegz')
    
    _cwd = os.getcwd()
    os.chdir(save_dir)

    logger.debug(f'Prepare to dock {docking_ligand} on {pdbid}')
    # 如果已有对接成功文件 跳过对接步骤
    if os.path.exists(dock_result_file) and not overwrite:     
        logger.debug('File %s is existed.' % dock_result_file)
        os.chdir(_cwd)                           
        return DockResultFile(dock_result_file, internal_ligand, internal_ligand_resnum, docking_ligand, precision)

    input_file = f'{pdbid}_{internal_ligand}_glide-dock_{docking_ligand}_{precision}.in'
    job_name = f'{pdbid}-{internal_ligand}-Glide-Dock-{docking_ligand}-{precision}'
    output_file = job_name + '_pv.maegz'
    glide_dock_config = [
        'FORCEFIELD %s\n' % FORCEFIELD,
        'GRIDFILE %s\n' % grid_file.file_path,
        'LIGANDFILE %s\n' % lig_file.file_path,
        'PRECISION %s\n' % precision
    ]
    if ligand_only:
        glide_dock_config.append('POSE_OUTTYPE ligandlib\n')
    if calc_rmsd:
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
        #logger.debug(f'{pdbid}-{docking_ligand} Docking Failed')
        raise RuntimeError(f'{pdbid}-{docking_ligand} Docking Failed')

    logger.debug(f'{pdbid}-{docking_ligand} Glide Docking Completed')
    logger.debug(f'Docking Result File: {dock_result_file} Saved.')
    os.chdir(_cwd)

    return DockResultFile(dock_result_file, internal_ligand, internal_ligand_resnum, docking_ligand, precision)

def calc_mmgbsa(maestrofile:DockResultFile, overwrite:bool=False) -> ComplexFile:
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
    
    prefix = maestrofile.file_prefix
    logger.debug('Prepare to Calculate MM-GB/SA Binding Energy : %s' % prefix)
    mmgbsa_result_file = prefix + '_mmgbsa.maegz'
    job_name = f'{prefix}_mmgbsa'

    # 已计算则跳过计算过程
    if os.path.exists(mmgbsa_result_file) and not overwrite:
        logger.debug('File %s is existed.' % mmgbsa_result_file)
        return ComplexFile(mmgbsa_result_file)

    launch(f'prime_mmgbsa -j {job_name} {maestrofile.file_path}')
    
    try:
        os.rename(f'{job_name}-out.maegz', mmgbsa_result_file)
    except Exception:
        logger.debug(f'{prefix} Prime MM-GB/SA Calculating Failed')
        return None

    logger.debug(f'MM-GB/SA Calculating Result File: {mmgbsa_result_file} Saved.')

    return ComplexFile(mmgbsa_result_file)

def calc_volume(recep_file:ReceptorFile, lig_file:LigandFile, overwrite:bool=False) -> ComplexFile:
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
        logger.debug('File %s is existed.' % sitemap_result_file)
        return ComplexFile(sitemap_result_file)

    launch(f'sitemap -sitebox 6 -keeplogs yes -ligmae {lig_file.file_path} -prot {recep_file.file_path} -j {pdbid}_sitemap')
    
    if not os.path.exists(sitemap_result_file):
        logger.debug(f'{pdbid} Sitemap Calculating Failed')
        return None

    logger.debug('Sitemap Calculating File: %s Saved.' % sitemap_result_file)
    return ComplexFile(sitemap_result_file)

def calc_admet(lig_file:LigandFile, overwrite:bool=False) -> LigandFile:
    '''
    计算化合物/配体的ADMET特征描述符

    Parameters  
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
        logger.debug('File %s is existed.' % admet_result_file)
        return LigandFile(admet_result_file)

    launch(f'qikprop -outname {prefix} {lig_file.file_path}')
    try:
        os.rename(prefix + '-out.mae', admet_result_file)
    except Exception:
        logger.debug(f'{prefix} ADMET Calculating Failed')
        return None

    logger.debug('Qikprop Calculation File: %s Saved.' % admet_result_file)
    return LigandFile(admet_result_file)
    