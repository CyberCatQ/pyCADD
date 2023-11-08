import os
import re
import logging
from time import sleep
from typing import Union

from pyCADD.Dynamic.template import LeapInput
from pyCADD.Density.base import Gauss
from pyCADD.utils.common import BaseFile
from pyCADD.utils.tool import _get_progress, makedirs_from_list, timeit

CWD = os.getcwd()
CPU_NUM = os.cpu_count() or 2
PMEMD = 'pmemd.cuda'
SANDER = f'mpirun -np {str(int(CPU_NUM / 2))} sander.MPI'

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
TEMPLATE_DIR = os.path.join(SCRIPT_DIR, 'template')

logger = logging.getLogger(__name__)


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
        raise RuntimeError(f'System call failed: {cmd}') from None


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


def _get_atom_lines(mol2file):
    '''
    从mol2文件中获取原子行信息

    Parameters
    ----------
    mol2file : str
        mol2文件路径

    Returns
    -------
    list
        原子行信息列表 
        [
            [原子序号, 原子名, x, y, z, 原子类型, 残基序号, 残基名, 电荷量],
            ...
        ]
    '''
    with open(mol2file) as f:
        lines = f.read().splitlines()
    start_idx = lines.index('@<TRIPOS>ATOM')
    end_idx = lines.index('@<TRIPOS>BOND')
    atom_lines = lines[start_idx+1:end_idx]
    return [line.split() for line in atom_lines]


def _merge_charge(charge_mol2, origin_mol2, output_mol2):
    '''
    合并mol2原始坐标与计算得到的电荷量

    Parameters
    ----------
    charge_mol2 : str
        计算得到电荷量的mol2文件路径(原子顺序需与origin_mol2一致)
    origin_mol2 : str
        原始mol2文件路径
    output_mol2 : str
        合并后的mol2文件路径

    '''
    mol2_line_fmt = "{:>7} {:<10} {:>10} {:>10} {:>10} {:<5} {:>4} {:<7} {:>10}\n"
    ori_atom_lines = _get_atom_lines(origin_mol2)
    chg_atom_lines = _get_atom_lines(charge_mol2)
    result_lines = ori_atom_lines.copy()
    for i, line in enumerate(result_lines):
        result_lines[i][8] = chg_atom_lines[i][8]
    result_lines = [mol2_line_fmt.format(*line) for line in result_lines]

    with open(origin_mol2) as f:
        lines = f.readlines()
    atom_start_idx = lines.index('@<TRIPOS>ATOM\n')
    atom_end_idx = lines.index('@<TRIPOS>BOND\n')

    with open(output_mol2, 'w') as f:
        f.writelines(lines[:atom_start_idx+1])
        f.writelines(result_lines)
        f.writelines(lines[atom_end_idx:])


def protein_prepare(protein_file: BaseFile, save_dir: Union[str, None] = None, keep_water: bool = False) -> BaseFile:
    '''
    预处理蛋白质PDB文件
        PDB文件格式化 for Amber | 去除原生H原子 | 使用rudece添加H原子 | 再次格式化

    Parameters
    ----------
    protein_file : BaseFile
        蛋白质PDB文件
    save_dir : str, optional
        保存路径 过程及结果文件保存至该目录 如为None则保存至当前目录
    keep_water : bool, optional
        是否保留输入结构中的水分子, 默认为 False

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

    if keep_water:
        _system_call(f'pdb4amber -i {file_path} -o tmp.pdb')
        _system_call(f'cat tmp.pdb | grep "HOH" > water.pdb')
        _system_call(
            f'cat {leap_file_path} | grep -v "END" > tmp.pdb && \
            cat water.pdb | sed "s/HOH/WAT/" >> tmp.pdb && \
            echo END >> tmp.pdb')
        _system_call(f'pdb4amber -i tmp.pdb -o {leap_file_path}')
        _system_call(f'rm tmp.pdb water.pdb')
    return BaseFile(leap_file_path)


def molecule_prepare_resp2(
        ligand_file: BaseFile,
        cpu_num: Union[int, None] = None,
        charge: Union[int, None] = None,
        multiplicity: Union[int, None] = None,
        solvent: Union[str, None] = None,
        save_dir: Union[str, None] = None,
        overwrite: bool = False,
        keep_origin_cood: bool = False) -> tuple:
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
    overwrite : bool, optional
        是否覆盖已存在文件 默认为False
    keep_origin_cood : bool, optional
        是否在输出结构中保留原始坐标 而不使用高斯结构优化的坐标 默认为False

    Returns
    -------
    BaseFile, BaseFile
        已计算电荷的mol2文件, frcmod文件
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

    # 参数生成过程文件
    pqr_file_path = os.path.join(save_dir, file_prefix + '.pqr')
    mol2_file_path = os.path.join(save_dir, file_prefix + '.mol2')
    frcmod_file_path = os.path.join(save_dir, file_prefix + '.frcmod')

    if os.path.exists(mol2_file_path) and not overwrite:
        _system_call(
            f'parmchk2 -i {mol2_file_path} -f mol2 -o {frcmod_file_path}')
        os.chdir(CWD)
        return BaseFile(mol2_file_path), BaseFile(frcmod_file_path)

    Gauss(file_path).set_system(cpu_num, '16GB')

    if os.popen('which obabel').read() == '':
        raise RuntimeError('Openbabel may not be installed.')
    print('Optimizing ligand structure...')
    # 高斯坐标优化与RESP2电荷计算
    # 完成时 生成.pqr文件
    _system_call(
        f'chmod 777 {script_resp2_path} && {script_resp2_path} {file_path} {charge} {multiplicity} {solvent}'
    )

    # Openbabel格式转换 pqr -> mol2 电荷信息传递至mol2文件
    _system_call(f'obabel -ipqr {pqr_file_path} -omol2 -O tmp.mol2')
    # Openbabel生成的mol2文件无法被parmchk2直接使用 需要antechamber转换
    # mol2文件带有RESP2电荷信息 但坐标已经高斯优化而改变
    _system_call(
        f'antechamber -fi mol2 -i tmp.mol2 -fo mol2 -o {mol2_file_path} && rm tmp.mol2')

    if keep_origin_cood:
        # 维持对接结果的坐标 并将RESP2电荷信息合并至mol2文件
        origin_file = ligand_file.file_path
        file_ext = ligand_file.file_ext

        if file_ext == 'pdb':
            _system_call(f'pdb4amber -i {file_path} -o origin.pdb')
            _system_call(
                f'antechamber -fi pdb -i origin.pdb -fo mol2 -o origin.mol2 && rm origin.pdb')
        else:
            _system_call(
                f'antechamber -fi {file_ext} -i {file_path} -fo mol2 -o origin.mol2')
        _merge_charge(mol2_file_path, 'origin.mol2', mol2_file_path)

    # 生成Amber Gaff Force Filed Parameters文件
    _system_call(
        f'parmchk2 -i {mol2_file_path} -f mol2 -o {frcmod_file_path}')

    os.chdir(CWD)

    return BaseFile(mol2_file_path), BaseFile(frcmod_file_path)


def molecule_prepare_bcc(
        ligand_file: BaseFile,
        charge: int,
        save_dir: Union[str, None] = None,
        overwrite: bool = False) -> tuple:
    '''
    使用AM1-BCC快速计算原子电荷 并完成配体预处理

    Parameters
    ----------
    ligand_file : BaseFile
        小分子文件
    charge : int
        电荷数
    save_dir : str, optional
        保存路径 过程及结果文件保存至该目录 如为None则保存至当前目录
    overwrite : bool, optional
        是否覆盖已存在的mol2结果文件 默认为False
    '''

    save_dir = save_dir if save_dir is not None else CWD
    file_path = ligand_file.file_path
    file_prefix = ligand_file.file_prefix

    # 参数生成过程文件
    mol2_file_path = os.path.join(save_dir, file_prefix + '.mol2')
    frcmod_file_path = os.path.join(save_dir, file_prefix + '.frcmod')

    if os.path.exists(mol2_file_path) and not overwrite:
        _system_call(
            f'parmchk2 -i {mol2_file_path} -f mol2 -o {frcmod_file_path}')
        return BaseFile(mol2_file_path), BaseFile(frcmod_file_path)

    _system_call(
        f'pdb4amber -i {ligand_file.file_path} -o tmp.pdb'
    )

    _system_call(
        f'antechamber -fi pdb -i tmp.pdb -fo mol2 -o {mol2_file_path} -c bcc -nc {charge} && rm tmp.pdb'
    )

    _system_call(
        f'parmchk2 -i {mol2_file_path} -f mol2 -o {frcmod_file_path}'
    )

    return BaseFile(mol2_file_path), BaseFile(frcmod_file_path)


def _creat_leap_inputfile(
        prefix: str,
        protein_file_path: str,
        ligand_file_path: Union[str, None] = None,
        frcmod_file_path: Union[str, None] = None,
        box_size: float = 12.0,
        save_dir: Union[str, None] = None) -> BaseFile:
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

    # tamplete_file_path = os.path.join(TEMPLATE_DIR, 'leap_template.in')
    input_file_path = os.path.join(save_dir, prefix + '_leap.in')

    input_file = LeapInput(
        protein_file_path=protein_file_path,
        ligand_file_path=ligand_file_path,
        frcmod_file_path=frcmod_file_path,
        file_prefix=prefix,
        box_size=box_size
    )
    input_file.save(input_file_path)
    return BaseFile(input_file_path)


def leap_prepare(prefix: str, ligand_file: BaseFile, frcmod_file: BaseFile, protein_file: BaseFile, box_size: float = 12.0, save_dir: Union[str, None] = None) -> None:
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
    box_size : float, optional
        水箱大小 默认为12.0 Angstrom
    save_dir : str, optional
        保存路径 默认为当前目录
    '''
    ligand_file_path = ligand_file.file_path
    frcmod_file_path = frcmod_file.file_path
    protein_file_path = protein_file.file_path

    save_dir = save_dir if save_dir is not None else CWD
    os.chdir(save_dir)

    leap_inputfile = _creat_leap_inputfile(
        prefix, ligand_file_path=ligand_file_path, frcmod_file_path=frcmod_file_path, protein_file_path=protein_file_path, box_size=box_size, save_dir=save_dir)
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


def leap_prepare_for_apo(prefix: str, protein_file: BaseFile, box_size: float = 12.0, save_dir: Union[str, None] = None) -> None:
    '''
    创建Apo晶体的LEaP输入文件并执行tleap命令

    Parameters
    -----------
    prefix : str
        生成文件的文件名前缀
    protein_file : BaseFile
        蛋白质PDB文件
    box_size : float, optional
        水箱大小 默认为12.0 Angstrom
    save_dir : str, optional
        保存路径 默认为当前目录
    '''

    protein_file_path = protein_file.file_path

    save_dir = save_dir if save_dir is not None else CWD
    os.chdir(save_dir)

    leap_inputfile = _creat_leap_inputfile(
        prefix=prefix, protein_file_path=protein_file_path, box_size=box_size, save_dir=save_dir)
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


def _trace_progress(output_file_path: str, step: int):
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


def _get_input_config(input_file: str) -> list:
    '''
    获取分子动力学模拟过程的输入文件配置

    Parameters
    -----------
    input_file : str
        分子动力学模拟过程的输入文件路径

    Returns
    -----------
    dict
        分子动力学模拟过程的输入文件配置
    '''
    output_list = []
    type_pattern = r'&(.*?)\n'
    config_pattern = "(?<=&{_type}\n).*"

    with open(input_file, 'r') as f:
        lines = f.readlines()

    config_str = "".join([line for line in lines])
    config_list = config_str.split('/')

    for index, config in enumerate(config_list):
        config = config.strip()
        if config == '' or config == 'END':
            continue
        _type = re.findall(type_pattern, config)[0]
        _config_list = re.findall(config_pattern.format(
            _type=_type), config, re.S)[0].split(',')
        config_dict = {"_index": index, "_type": _type}
        config_dict.update(
            {k.strip(): v.strip()
             for k, v in [item.split('=') for item in _config_list]}
        )
        output_list.append(config_dict)
    return output_list


class BaseProcess:
    def __init__(self, input_file: BaseFile, process_name: str, **kwargs) -> None:
        self.input_file = input_file
        self.name = self.process_name = process_name
        self.cfg = _get_input_config(input_file.file_path)
        for k, v in kwargs.items():
            setattr(self, k, v)

    def run(self, **kwargs) -> None:
        raise NotImplementedError


class MDProcess(BaseProcess):
    def __init__(self, input_file: BaseFile, process_name: str, **kwargs) -> None:

        super().__init__(input_file, process_name, **kwargs)
        self.control_cfg = [
            cfg for cfg in self.cfg if cfg['_type'] == 'cntrl'][0]
        self.is_minimize = bool(int(self.control_cfg.get('imin') or 0))
        self.is_restrained = bool(int(self.control_cfg.get('ntr') or 0))
        self.is_production = False
        self.total_step = int(self.control_cfg.get('nstlim') or self.control_cfg.get('maxcyc'))
        self.step_size = float(self.control_cfg.get('dt')) if not self.is_minimize else 1

        self.toplogy_file = None
        self.inpcrd_file = None
        self.mdout_file_path = None
        self.mdcrd_file_path = None
        self.mdrst_file_path = None

    def run(self, toplogy_file: BaseFile, inpcrd_file: BaseFile, reference_file: Union[BaseFile, None] = None, save_dir: Union[str, None] = None, nohup: bool = False) -> 'MDProcess':
        '''
        Run MD process.

        Parameters
        -----------
        save_dir : str, optional
            Save directory. Default is CWD.
        '''
        self.toplogy_file = toplogy_file
        self.inpcrd_file = inpcrd_file
        save_dir = save_dir if save_dir is not None else CWD
        os.makedirs(save_dir, exist_ok=True)

        self.mdout_file_path = os.path.join(
            save_dir, self.process_name + '.out')
        self.mdcrd_file_path = os.path.join(
            save_dir, self.process_name + '.nc')
        self.mdrst_file_path = os.path.join(
            save_dir, self.process_name + '.rst')
        simulation_time = self.total_step * self.step_size / 1000

        if not self.is_minimize:
            logger.info(f'Simulation total steps: {self.total_step}')
            logger.info(f'Simulation step size: {self.step_size} ps')
            logger.info(f'Simulation total time: {simulation_time} ns')

        self.cmd = f"{PMEMD} -O"
        self.cmd += f" -i {self.input_file.file_path}"
        self.cmd += f" -o {self.mdout_file_path}"
        self.cmd += f" -c {self.inpcrd_file.file_path}"
        self.cmd += f" -p {self.toplogy_file.file_path}"
        self.cmd += f" -r {self.mdrst_file_path}"

        if not self.is_minimize:
            self.cmd += f" -x {self.mdcrd_file_path}"
        else:
            if reference_file is not None:
                self.cmd += f" -ref {reference_file.file_path}"

        if nohup:
            self.cmd = f'nohup {self.cmd} > /dev/null 2>&1 &'

        logger.debug(
            f'Process {self.process_name} running command: {self.cmd}')
        logger.info(f'Running process {self.process_name} ...')
        _system_call(self.cmd)
        if nohup:
            _trace_progress(self.mdout_file_path, self.total_step)
        # try:
        #     _system_call(self.cmd)
        #     if nohup:
        #         _trace_progress(self.mdout_file_path, self.total_step)
        # except RuntimeError:
        #     logger.warning(
        #     'GPU version of PMEMD failed, trying CPU version sander.MPI ...')
        #     cmd_cpu = self.cmd.replace(PMEMD, SANDER)
        #     logger.info(f'Using CMD: {cmd_cpu}')
        #     _system_call(cmd_cpu)
        #     if not self.is_minimize:
        #         _trace_progress(self.mdout_file_path, self.total_step)

        logger.info(f"Process {self.process_name} finished.")

        return self


class MinimizeProcess(MDProcess):
    def __init__(self, input_file: BaseFile, process_name: str, **kwargs) -> None:
        super().__init__(input_file, process_name, **kwargs)


class NPTProcess(MDProcess):
    def __init__(self, input_file: BaseFile, process_name: str, is_production: bool = False, **kwargs) -> None:
        super().__init__(input_file, process_name, **kwargs)
        self.is_production = is_production


class NVTProcess(MDProcess):
    def __init__(self, input_file: BaseFile, process_name: str, is_production: bool = False, **kwargs) -> None:
        super().__init__(input_file, process_name, **kwargs)
        self.is_production = is_production


@timeit
def _run_simulation(
        comsolvate_topfile: BaseFile,
        comsolvate_crdfile: BaseFile,
        process_list: list,
        save_dir: Union[str, None] = None) -> None:
    '''
    执行分子动力学模拟工作流
    '''
    save_dir = save_dir if save_dir is not None else CWD

    for index, process in enumerate(process_list):

        process_output_dir = os.path.join(save_dir, process.process_name)
        inpcrd_file = comsolvate_crdfile if index == 0 else BaseFile(
            finished_process.mdrst_file_path)
        reference_file = inpcrd_file if process.is_restrained else None
        nohup = True if process.is_production else False
        finished_process = process.run(
            toplogy_file=comsolvate_topfile,
            inpcrd_file=inpcrd_file,
            save_dir=process_output_dir,
            reference_file=reference_file,
            nohup=nohup
        )
        sleep(1)

    # step_a_dir = os.path.join(save_dir, 'step_a')
    # step_b_dir = os.path.join(save_dir, 'step_b')
    # step_c_dir = os.path.join(save_dir, 'step_c')
    # step_nvt_dir = os.path.join(save_dir, 'step_nvt')
    # step_npt_dir = os.path.join(save_dir, 'step_npt')

    # nvt_cfg = _get_input_config(step_nvt_inputfile.file_path)
    # npt_cfg = _get_input_config(step_npt_inputfile.file_path)

    # step_num = int(npt_cfg['nstlim'])
    # step_size = float(npt_cfg['dt'])
    # simulate_time = step_num * step_size / 1000

    # logger.info(f'Simulate steps: {step_num}')
    # logger.info(f'Simulate step size: {step_size}')
    # logger.info(f'Simulate time: {simulate_time} ns')

    # makedirs_from_list(
    #     [
    #         step_a_dir,
    #         step_b_dir,
    #         step_c_dir,
    #         step_nvt_dir,
    #         step_npt_dir
    #     ]
    # )

    # step_a_outfile = os.path.join(step_a_dir, 'step_a.out')
    # step_b_outfile = os.path.join(step_b_dir, 'step_b.out')
    # step_c_outfile = os.path.join(step_c_dir, 'step_c.out')
    # step_nvt_outfile = os.path.join(step_nvt_dir, 'step_nvt.out')
    # step_npt_outfile = os.path.join(step_npt_dir, 'step_npt.out')

    # step_a_rstfile = os.path.join(step_a_dir, 'step_a.rst')
    # step_b_rstfile = os.path.join(step_b_dir, 'step_b.rst')
    # step_c_rstfile = os.path.join(step_c_dir, 'step_c.rst')
    # step_nvt_rstfile = os.path.join(step_nvt_dir, 'step_nvt.rst')
    # step_npt_rstfile = os.path.join(step_npt_dir, 'step_npt.rst')

    # step_nvt_crdfile = os.path.join(step_nvt_dir, 'step_nvt.crd')
    # step_npt_crdfile = os.path.join(step_npt_dir, 'step_npt.crd')

    # # Get more information from AMBER USER PROFILE
    # step_a_cmd = f'{PMEMD} -O '
    # step_a_cmd += f'-i {step_a_inputfile.file_path} '
    # step_a_cmd += f'-o {step_a_outfile} '
    # step_a_cmd += f'-c {comsolvate_crdfile.file_path} -p {comsolvate_topfile.file_path} '
    # step_a_cmd += f'-r {step_a_rstfile} '
    # step_a_cmd += f'-ref {comsolvate_crdfile.file_path}'

    # step_b_cmd = f'{PMEMD} -O '
    # step_b_cmd += f'-i {step_b_inputfile.file_path} '
    # step_b_cmd += f'-o {step_b_outfile} '
    # step_b_cmd += f'-c {step_a_rstfile} -p {comsolvate_topfile.file_path} '
    # step_b_cmd += f'-r {step_b_rstfile} '
    # step_b_cmd += f'-ref {step_a_rstfile}'

    # step_c_cmd = f'{PMEMD} -O '
    # step_c_cmd += f'-i {step_c_inputfile.file_path} '
    # step_c_cmd += f'-o {step_c_outfile} '
    # step_c_cmd += f'-c {step_b_rstfile} -p {comsolvate_topfile.file_path} '
    # step_c_cmd += f'-r {step_c_rstfile} '

    # step_nvt_cmd = f'{PMEMD} -O '
    # step_nvt_cmd += f'-i {step_nvt_inputfile.file_path} '
    # step_nvt_cmd += f'-o {step_nvt_outfile} '
    # step_nvt_cmd += f'-c {step_c_rstfile} -p {comsolvate_topfile.file_path} '
    # step_nvt_cmd += f'-r {step_nvt_rstfile} '
    # step_nvt_cmd += f'-x {step_nvt_crdfile}'

    # step_nvt_cmd_cpu = f'{SANDER} -O '
    # step_nvt_cmd_cpu += f'-i {step_nvt_inputfile.file_path} '
    # step_nvt_cmd_cpu += f'-o {step_nvt_outfile} '
    # step_nvt_cmd_cpu += f'-c {step_c_rstfile} -p {comsolvate_topfile.file_path} '
    # step_nvt_cmd_cpu += f'-r {step_nvt_rstfile} '
    # step_nvt_cmd_cpu += f'-x {step_nvt_crdfile}'

    # step_npt_cmd = f'nohup {PMEMD} -O '
    # step_npt_cmd += f'-i {step_npt_inputfile.file_path} '
    # step_npt_cmd += f'-o {step_npt_outfile} '
    # step_npt_cmd += f'-c {step_nvt_rstfile} -p {comsolvate_topfile.file_path} '
    # step_npt_cmd += f'-r {step_npt_rstfile} '
    # step_npt_cmd += f'-x {step_npt_crdfile} > /dev/null 2>&1 &'

    # logger.info('Running Minimize Progress A...')
    # _system_call(step_a_cmd)
    # logger.info('Running Minimize Progress B...')
    # _system_call(step_b_cmd)
    # logger.info('Running Minimize Progress C...')
    # _system_call(step_c_cmd)
    # logger.info('Running 100ps Heating Progress...')

    # try:
    #     _system_call(step_nvt_cmd)
    # except RuntimeError:
    #     logger.warning(
    #         'GPU version of PMEMD failed, trying CPU version sander.MPI ...')
    #     logger.info(f'Using CMD: {SANDER}')
    #     _system_call(step_nvt_cmd_cpu)
    #     _trace_progress(step_nvt_outfile, 50000)

    # logger.info(f'Running {simulate_time} ns Molecular Dynamics Simulation...')
    # _system_call(step_npt_cmd)
    # _trace_progress(step_npt_outfile, step_num)

    # logger.info('Simulation Finished.')


def _creat_energy_inputfile(job_type: str, startframe: int, endframe: int, interval: int, decomp: bool = False, save_dir: Union[str, None] = None) -> BaseFile:
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


@timeit
def _run_energy_calculation(
    input_file: BaseFile, comsolvate_topfile: BaseFile, com_topfile: BaseFile,
    receptor_topfile: BaseFile, ligand_topfile: BaseFile, traj_file: BaseFile,
    output_filepath: Union[str, None] = None, decom_output_filepath: Union[str, None] = None, cpu_num: Union[int, None] = None
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
