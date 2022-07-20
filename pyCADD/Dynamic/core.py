import os
import re
from time import sleep

from pyCADD.Gauss.base import Gauss
from pyCADD.utils.common import BaseFile
from pyCADD.utils.tool import _get_progress, makedirs_from_list

CWD = os.getcwd()
CPU_NUM = os.cpu_count()
SANDER = 'pmemd.cuda'

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
TEMPLATE_DIR = os.path.join(SCRIPT_DIR, 'template')


def _system_call(cmd: str) -> None:
    '''
    系统命令调用

    Parameters
    ----------
    cmd : str
        命令行

    '''

    return_code = os.system(cmd)
    if return_code != 0:
        raise RuntimeError(f'System call failed: {cmd}')


def _convert_mae_to_pdb(mae_file_path: str) -> None:
    '''
    转换MAE文件为PDB文件

    Parameters
    ----------
    mae_file_path : str
        mae文件路径
    '''
    script_path = os.path.join(SCRIPT_DIR, 'mae2pdb.sh')
    _system_call(f'{script_path} {mae_file_path}')


def _creat_file_from_template(template_file_path: str, output_file_path: str, **kwarg):
    '''
    从模板文件中读取模板信息
    并替换为指定内容 从而创建定制文件
    '''
    with open(template_file_path, 'r') as f:
        template = f.read()

    output_file_text = template.format(**kwarg)

    with open(output_file_path, 'w') as f:
        f.write(output_file_text)

    return BaseFile(output_file_path)


def protein_prepare(protein_file: BaseFile, save_dir: str = None) -> BaseFile:
    '''
    预处理蛋白质PDB文件
        PDB文件格式化 for Amber | 去除原生H原子 | 使用rudece添加H原子 | 再次格式化

    Parameters
    ----------
    protein_file : BaseFile
        蛋白质PDB文件
    save_dir : str, optional
        保存路径 过程及结果文件保存至该目录 如为None则保存至当前目录

    Returns
    -------
    str
        处理完成的蛋白质PDB文件
    '''

    save_dir = save_dir if save_dir is not None else CWD
    file_path = protein_file.file_path
    file_prefix = protein_file.file_prefix

    dry_file_path = os.path.join(save_dir, file_prefix + '_dry.pdb')
    noH_file_path = os.path.join(save_dir, file_prefix + '_noH.pdb')
    leap_file_path = os.path.join(save_dir, file_prefix + '_leap.pdb')

    # 去除水分子
    _system_call(
        f'pdb4amber -i {file_path} -o {dry_file_path} -p --dry --add-missing-atoms ')
    # 去除原生H原子
    _system_call(f'pdb4amber -i {dry_file_path} -o {noH_file_path} -y')
    # 重新格式化
    _system_call(f'pdb4amber -i {noH_file_path} -o {leap_file_path}')

    return BaseFile(leap_file_path)


def molecule_prepare(
        ligand_file: BaseFile,
        cpu_num: int = None,
        charge: int = None,
        multiplicity: int = None,
        solvent: str = None,
        save_dir: str = None) -> tuple:
    '''
    预处理小分子文件
        高斯坐标优化与RESP2(0.5)电荷计算

    Parameters
    ----------
    ligand_file : BaseFile
        小分子文件
    cpu_num : int, optional
        CPU数量 默认为CPU总数
    charge : int, optional
        电荷数 默认为电荷数为0
    multiplicity : int, optional
        自旋多重度 默认为复数为1
    solvent : str, optional
        计算溶剂 默认为水
    save_dir : str, optional
        保存路径 过程及结果文件保存至该目录 如为None则保存至当前目录
    '''

    save_dir = save_dir if save_dir is not None else CWD
    file_path = ligand_file.file_path
    file_prefix = ligand_file.file_prefix

    os.chdir(save_dir)

    cpu_num = cpu_num if cpu_num is not None else CPU_NUM
    charge = charge if charge is not None else 0
    multiplicity = multiplicity if multiplicity is not None else 1
    solvent = solvent if solvent is not None else 'water'
    script_resp2_path = os.path.join(SCRIPT_DIR, 'RESP2.sh')

    Gauss(file_path).set_system(cpu_num, '16GB')

    print('Optimizing ligand structure...')
    # 高斯坐标优化与RESP2电荷计算
    # 完成时 生成.pqr文件
    _system_call(
        f'chmod 777 {script_resp2_path} && {script_resp2_path} {file_path} {charge} {multiplicity} {solvent}'
    )

    # 参数生成过程文件
    pqr_file_path = os.path.join(save_dir, file_prefix + '.pqr')
    mol2_file_path = os.path.join(save_dir, file_prefix + '.mol2')

    # Openbabel格式转换 pqr -> mol2 电荷信息传递至mol2文件
    _system_call(f'obabel -ipqr {pqr_file_path} -omol2 -O tmp.mol2')
    # Openbabel生成的mol2文件无法被parmchk2直接使用 需要antechamber转换
    # mol2文件带有RESP2电荷信息
    _system_call(
        f'antechamber -fi mol2 -i tmp.mol2 -fo mol2 -o {mol2_file_path} && rm tmp.mol2')

    frcmod_file_path = os.path.join(save_dir, file_prefix + '.frcmod')

    # 生成Amber Gaff Force Filed Parameters文件
    _system_call(
        f'parmchk2 -i {mol2_file_path} -f mol2 -o {frcmod_file_path}')

    os.chdir(CWD)

    return BaseFile(mol2_file_path), BaseFile(frcmod_file_path)


def _creat_leap_inputfile(
        prefix: str,
        ligand_file_path: str,
        frcmod_file_path: str,
        protein_file_path: str,
        save_dir: str = None) -> BaseFile:
    '''
    创建LEaP输入文件

    Parameters
    -----------
    prefix : str
        文件名前缀
    ligand_file_path : str
        小分子文件路径
    frcmod_file_path : str
        Amber Parameters文件路径
    protein_file_path : str
        蛋白质PDB文件路径
    save_dir : str, optional
        保存路径
    '''
    save_dir = save_dir if save_dir is not None else CWD

    tamplete_file_path = os.path.join(TEMPLATE_DIR, 'leap_template.in')
    input_file_path = os.path.join(save_dir, prefix + '_leap.in')

    input_file = _creat_file_from_template(
        tamplete_file_path, input_file_path,
        ligand_file_path=ligand_file_path,
        frcmod_file_path=frcmod_file_path,
        protein_file_path=protein_file_path,
        pro_lig='{pro lig}',
        file_prefix=prefix
    )

    return input_file


def leap_prepare(prefix: str, ligand_file: BaseFile, frcmod_file: BaseFile, protein_file: BaseFile, save_dir: str = None) -> None:
    '''
    创建LEaP输入文件并执行tleap命令

    Parameters
    -----------
    prefix : str
        生成文件的文件名前缀
    ligand_file : BaseFile
        小分子文件
    frcmod_file : BaseFile
        Amber Parameters文件
    protein_file : BaseFile
        蛋白质PDB文件
    save_dir : str, optional
        保存路径 默认为当前目录
    '''
    ligand_file_path = ligand_file.file_path
    frcmod_file_path = frcmod_file.file_path
    protein_file_path = protein_file.file_path

    save_dir = save_dir if save_dir is not None else CWD
    os.chdir(save_dir)

    leap_inputfile = _creat_leap_inputfile(
        prefix, ligand_file_path, frcmod_file_path, protein_file_path, save_dir)
    # 调用tleap命令
    _system_call(f'tleap -f {leap_inputfile.file_path}')

    # 生成水箱复合物PDB文件
    comsolvate_topfile_path = os.path.join(
        save_dir, prefix + '_comsolvate.prmtop')
    comsolvate_crdfile_path = os.path.join(
        save_dir, prefix + '_comsolvate.inpcrd')
    comsolvate_pdbfile_path = os.path.join(
        save_dir, prefix + '_comsolvate.pdb')
    _system_call(
        f'ambpdb -p {comsolvate_topfile_path} < {comsolvate_crdfile_path} > {comsolvate_pdbfile_path}')

    os.chdir(CWD)


def _get_water_resnum(comsolvate_pdbfile: BaseFile) -> list:
    '''
    获取水箱复合物的全部水分子Residue Number

    Parameters
    -----------
    comsolvate_pdbfile : BaseFile
        水箱复合物PDB文件

    Returns
    -----------
    list
        水分子Residue Number列表
    '''
    return os.popen(
        "cat %s | grep WAT -w | awk '{print $5}'" % comsolvate_pdbfile.file_path).read().splitlines()


def _creat_md_inputfile(
        water_resnum: list,
        step_num: int = 50000000,
        step_length: float = 0.002,
        save_dir: str = None) -> tuple:
    '''
    创建分子动力学模拟过程的输入文件
        Step A: 约束主链能量最小化
        Step B: 约束非主链能量最小化
        Step C: 无约束能量最小化
        Step nvt: 体系恒容恒压加热100ps
        Step npt: 主模拟

    Parameters
    -----------
    water_resnum : list
        水分子Residue Number列表
    step_num : int
        模拟步数
    step_length : float
        模拟步长
    save_dir : str, optional
        保存路径 默认为当前目录

    Returns
    -----------
    tuple[BaseFile, BaseFile, BaseFile, BaseFile, BaseFile]
        分子动力学模拟过程的输入文件
        step A, B, C, nvt, npt
    '''
    save_dir = save_dir if save_dir is not None else CWD

    step_a_temfile_path = os.path.join(TEMPLATE_DIR, 'step_a.in')
    step_b_temfile_path = os.path.join(TEMPLATE_DIR, 'step_b.in')
    step_c_temfile_path = os.path.join(TEMPLATE_DIR, 'step_c.in')
    step_nvt_temfile_path = os.path.join(TEMPLATE_DIR, 'step_nvt.in')
    step_npt_temfile_path = os.path.join(TEMPLATE_DIR, 'step_npt.in')

    water_resnum_start = int(water_resnum[0])
    water_resnum_end = int(water_resnum[-1])

    # setp A input file
    step_a_inputfile_path = os.path.join(save_dir, 'step_a.in')
    # 限制主链能量最小化
    step_a_restraintmask = f"':1-{water_resnum_start-1}'"
    step_a_inputfile = _creat_file_from_template(
        step_a_temfile_path, step_a_inputfile_path,
        restraintmask=step_a_restraintmask
    )

    # step B input file
    step_b_inputfile_path = os.path.join(save_dir, 'step_b.in')
    # 限制非主链能量最小化
    step_b_restraintmask = f"':{water_resnum_start}-{water_resnum_end}'"
    step_b_inputfile = _creat_file_from_template(
        step_b_temfile_path, step_b_inputfile_path,
        restraintmask=step_b_restraintmask
    )

    # step C input file
    step_c_inputfile_path = os.path.join(save_dir, 'step_c.in')
    step_c_inputfile = _creat_file_from_template(
        step_c_temfile_path, step_c_inputfile_path
    )

    # step nvt input file
    step_nvt_inputfile_path = os.path.join(save_dir, 'step_nvt.in')
    step_nvt_inputfile = _creat_file_from_template(
        step_nvt_temfile_path, step_nvt_inputfile_path,
    )

    # step npt input file
    step_npt_inputfile_path = os.path.join(save_dir, 'step_npt.in')
    step_npt_inputfile = _creat_file_from_template(
        step_npt_temfile_path, step_npt_inputfile_path,
        step_num=step_num, step_length=step_length
    )

    return step_a_inputfile, step_b_inputfile, step_c_inputfile, step_nvt_inputfile, step_npt_inputfile


def _trace_progress(output_file_path: str, step: int = 50000000):
    '''
    追踪分子动力学模拟进程

    Parameters
    -----------
    output_file_path : str
        分子动力学模拟过程的输出文件路径
    step : int, optional
        模拟步数 默认为50000000
    '''
    progress, taskID = _get_progress(
        'Molecule Dynamics Simulation', 'bold cyan', total=step)
    progress.start()
    progress.start_task(taskID)
    while not progress.finished:
        _current = os.popen(
            f'tail -n 10 {output_file_path} | grep NSTEP').read()
        _finished = os.popen(
            f'tail -n 20 {output_file_path} | grep Final').read()
        if _current:
            current_step = re.findall(r'\d+', _current)[0]
            progress.update(taskID, completed=int(current_step))
        elif _finished:
            progress.update(taskID, completed=step)
            break
        sleep(1)
    progress.stop()


def _run_simulation(
        comsolvate_topfile: BaseFile,
        comsolvate_crdfile: BaseFile,
        step_a_inputfile: BaseFile,
        step_b_inputfile: BaseFile,
        step_c_inputfile: BaseFile,
        step_nvt_inputfile: BaseFile,
        step_npt_inputfile: BaseFile,
        save_dir: str = None) -> None:
    '''
    执行分子动力学模拟工作流
    '''
    save_dir = save_dir if save_dir is not None else CWD

    step_a_dir = os.path.join(save_dir, 'step_a')
    step_b_dir = os.path.join(save_dir, 'step_b')
    step_c_dir = os.path.join(save_dir, 'step_c')
    step_nvt_dir = os.path.join(save_dir, 'step_nvt')
    step_npt_dir = os.path.join(save_dir, 'step_npt')

    makedirs_from_list(
        [
            step_a_dir,
            step_b_dir,
            step_c_dir,
            step_nvt_dir,
            step_npt_dir
        ]
    )

    step_a_outfile = os.path.join(step_a_dir, 'step_a.out')
    step_b_outfile = os.path.join(step_b_dir, 'step_b.out')
    step_c_outfile = os.path.join(step_c_dir, 'step_c.out')
    step_nvt_outfile = os.path.join(step_nvt_dir, 'step_nvt.out')
    step_npt_outfile = os.path.join(step_npt_dir, 'step_npt.out')

    step_a_rstfile = os.path.join(step_a_dir, 'step_a.rst')
    step_b_rstfile = os.path.join(step_b_dir, 'step_b.rst')
    step_c_rstfile = os.path.join(step_c_dir, 'step_c.rst')
    step_nvt_rstfile = os.path.join(step_nvt_dir, 'step_nvt.rst')
    step_npt_rstfile = os.path.join(step_npt_dir, 'step_npt.rst')

    step_nvt_crdfile = os.path.join(step_nvt_dir, 'step_nvt.crd')
    step_npt_crdfile = os.path.join(step_npt_dir, 'step_npt.crd')

    # Get more information from AMBER USER PROFILE
    step_a_cmd = f'{SANDER} -O '
    step_a_cmd += f'-i {step_a_inputfile.file_path} '
    step_a_cmd += f'-o {step_a_outfile} '
    step_a_cmd += f'-c {comsolvate_crdfile.file_path} -p {comsolvate_topfile.file_path} '
    step_a_cmd += f'-r {step_a_rstfile} '
    step_a_cmd += f'-ref {comsolvate_crdfile.file_path}'

    step_b_cmd = f'{SANDER} -O '
    step_b_cmd += f'-i {step_b_inputfile.file_path} '
    step_b_cmd += f'-o {step_b_outfile} '
    step_b_cmd += f'-c {step_a_rstfile} -p {comsolvate_topfile.file_path} '
    step_b_cmd += f'-r {step_b_rstfile} '
    step_b_cmd += f'-ref {step_a_rstfile}'

    step_c_cmd = f'{SANDER} -O '
    step_c_cmd += f'-i {step_c_inputfile.file_path} '
    step_c_cmd += f'-o {step_c_outfile} '
    step_c_cmd += f'-c {step_b_rstfile} -p {comsolvate_topfile.file_path} '
    step_c_cmd += f'-r {step_c_rstfile} '

    step_nvt_cmd = f'{SANDER} -O '
    step_nvt_cmd += f'-i {step_nvt_inputfile.file_path} '
    step_nvt_cmd += f'-o {step_nvt_outfile} '
    step_nvt_cmd += f'-c {step_c_rstfile} -p {comsolvate_topfile.file_path} '
    step_nvt_cmd += f'-r {step_nvt_rstfile} '
    step_nvt_cmd += f'-x {step_nvt_crdfile}'

    step_npt_cmd = f'nohup {SANDER} -O '
    step_npt_cmd += f'-i {step_npt_inputfile.file_path} '
    step_npt_cmd += f'-o {step_npt_outfile} '
    step_npt_cmd += f'-c {step_nvt_rstfile} -p {comsolvate_topfile.file_path} '
    step_npt_cmd += f'-r {step_npt_rstfile} '
    step_npt_cmd += f'-x {step_npt_crdfile} > /dev/null 2>&1 &'

    print('Running Minimize Progress A...')
    _system_call(step_a_cmd)
    print('Running Minimize Progress B...')
    _system_call(step_b_cmd)
    print('Running Minimize Progress C...')
    _system_call(step_c_cmd)
    print('Running 100ps Heating Progress...')
    _system_call(step_nvt_cmd)
    print('Running Molecular Dynamics Simulation...')
    _system_call(step_npt_cmd)
    _trace_progress(step_npt_outfile)

    print('Simulation Finished.')


def _creat_energy_inputfile(job_type: str, startframe: int, endframe: int, interval: int, decomp:bool=False, save_dir: str = None) -> BaseFile:
    '''
    创建能量计算相关输入文件

    Parameters
    ----------
    job_type : str
        能量计算类型，可选值为：
        'pb/gb': 同时进行MM-PB/GBSA自由能计算
        'gb': 进行MM-GBSA自由能计算
        'nmode': 进行 Normal Mode 熵变计算
    startframe : int
        计算分析起始帧
    endframe : int
        计算分析结束帧
    interval : int
        计算分析帧间隔
    decomp : bool, optional
        是否进行能量分解计算 默认为False
    save_dir : str, optional
        计算结果保存目录 默认为当前目录

    Returns
    -------
    BaseFile
        能量计算输入文件
    '''
    save_dir = save_dir if save_dir is not None else CWD

    if job_type == 'pb/gb':
        if decomp:
            template_file = os.path.join(TEMPLATE_DIR, 'mmpb_gbsa_decom.in')
        else:
            template_file = os.path.join(TEMPLATE_DIR, 'mmpb_gbsa_only.in')
    elif job_type == 'gb':
        if decomp:
            template_file = os.path.join(TEMPLATE_DIR, 'mmgbsa_decom.in')
        else:
            template_file = os.path.join(TEMPLATE_DIR, 'mmgbsa_only.in')
    elif job_type == 'nmode':
        template_file = os.path.join(TEMPLATE_DIR, 'nmode.in')
    else:
        raise ValueError(f'Invalid job_type: {job_type}')

    output_file = os.path.join(save_dir, os.path.basename(template_file))
    input_file = _creat_file_from_template(
        template_file,
        output_file,
        startframe=startframe,
        endframe=endframe,
        interval=interval
    )

    return input_file


def _run_energy_calculation(
    input_file: BaseFile, comsolvate_topfile: BaseFile, com_topfile: BaseFile,
    receptor_topfile: BaseFile, ligand_topfile: BaseFile, traj_file: BaseFile,
    output_filepath: str = None, decom_output_filepath: str = None, cpu_num: int = None
    ) -> None:

    cpu_num = cpu_num if cpu_num is not None else CPU_NUM
    energy_cmd = f'mpirun -np {cpu_num} MMPBSA.py.MPI -O '
    energy_cmd += f'-i {input_file.file_path} '
    energy_cmd += f'-sp {comsolvate_topfile.file_path} '
    energy_cmd += f'-cp {com_topfile.file_path} '
    energy_cmd += f'-rp {receptor_topfile.file_path} '
    energy_cmd += f'-lp {ligand_topfile.file_path} '
    energy_cmd += f'-y {traj_file.file_path} '

    if output_filepath is not None:
        energy_cmd += f'-o {output_filepath} '
    if decom_output_filepath is not None:
        energy_cmd += f'-do {decom_output_filepath} '

    energy_cmd += f'> {input_file.file_prefix}.log '

    print('Running Energy Calculation...')
    _system_call(energy_cmd)
    print('Energy Calculation Finished.')
