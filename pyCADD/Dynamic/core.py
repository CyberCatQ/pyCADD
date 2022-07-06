import os

from pyCADD.utils.common import BaseFile
from pyCADD.Gauss.base import Gauss

CWD = os.getcwd()
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CPU_NUM = os.cpu_count()


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

    save_dir = save_dir if save_dir is None else CWD
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
    # 添加H原子并重新格式化
    _system_call(f'pdb4amber -i {noH_file_path} -o {leap_file_path}')

    return BaseFile(leap_file_path)


def ligand_prepare(
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

    save_dir = save_dir if save_dir is None else CWD
    file_path = ligand_file.file_path
    file_prefix = ligand_file.file_prefix

    os.system(f'cd {save_dir}')

    cpu_num = cpu_num if cpu_num is not None else CPU_NUM
    charge = charge if charge is not None else 0
    multiplicity = multiplicity if multiplicity is not None else 1
    solvent = solvent if solvent is not None else 'water'

    Gauss(file_path).set_system(cpu_num, '16GB')

    script_resp2_path = os.path.join(SCRIPT_DIR, 'RESP2.sh')

    print('Optimizing ligand structure...')
    # 高斯坐标优化与RESP2电荷计算
    _system_call(
        f'chmod 777 {script_resp2_path} && {script_resp2_path} {file_path} {charge} {multiplicity} {solvent}'
    )

    pqr_file_path = os.path.join(save_dir, file_prefix + '.pqr')
    mol2_file_path = os.path.join(save_dir, file_prefix + '.mol2')
    prepin_file_path = os.path.join(save_dir, file_prefix + '.prepin')
    pdb_file_path = os.path.join(save_dir, file_prefix + '_out.pdb')
    frcmod_file_path = os.path.join(save_dir, file_prefix + '.frcmod')

    # 优化及计算完成的结构转换为mol2格式
    _system_call(f'obabel -ipqr {pqr_file_path} -omol2 -O {mol2_file_path}')

    # 生成Amber pripi文件
    _system_call(
        f'antechamber -fi mol2 -i {mol2_file_path} -fo prepi -o {prepin_file_path}')

    # 生成Amber pdb文件
    _system_call(
        f'antechamber -fi mol2 -i {mol2_file_path} -o {pdb_file_path} -fo pdb')

    # 生成Amber Parameters文件
    _system_call(
        f'parmchk2 -i {prepin_file_path} -f prepi -o {frcmod_file_path}')

    os.system(f'cd {CWD}')

    return BaseFile(pdb_file_path), BaseFile(prepin_file_path), BaseFile(frcmod_file_path)


def _creat_leap_inputfile(
        prefix: str,
        ligand_file_path: str,
        prepin_file_path: str,
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
    prepin_file_path : str
        Amber prepi文件路径
    frcmod_file_path : str
        Amber Parameters文件路径
    protein_file_path : str
        蛋白质PDB文件路径
    save_dir : str, optional
        保存路径
    '''
    save_dir = save_dir if save_dir is None else CWD

    tamplete_file_path = os.path.join(SCRIPT_DIR, 'leap_template.in')
    with open(tamplete_file_path, 'r') as f:
        template = f.read()

    input_file_text = template.format(
        ligand_file_path=ligand_file_path,
        prepin_file_path=prepin_file_path,
        frcmod_file_path=frcmod_file_path,
        protein_file_path=protein_file_path,
        pro_lig='{pro lig}',
        file_prefix=prefix
    )

    input_file_path = os.path.join(save_dir, prefix + '_leap.in')
    with open(input_file_path, 'w') as f:
        f.write(input_file_text)

    return BaseFile(input_file_path)


def leap_prepare(prefix: str, ligand_file: BaseFile, prepin_file: BaseFile, frcmod_file: BaseFile, protein_file: BaseFile, save_dir: str = None) -> None:
    '''
    创建LEaP输入文件并执行tleap命令

    Parameters
    -----------
    prefix : str
        生成文件的文件名前缀
    ligand_file : BaseFile
        小分子文件
    prepin_file : BaseFile
        Amber prepi文件
    frcmod_file : BaseFile
        Amber Parameters文件
    protein_file : BaseFile
        蛋白质PDB文件
    save_dir : str, optional
        保存路径 默认为当前目录
    '''
    ligand_file_path = ligand_file.file_path
    prepin_file_path = prepin_file.file_path
    frcmod_file_path = frcmod_file.file_path
    protein_file_path = protein_file.file_path

    save_dir = save_dir if save_dir is None else CWD
    os.system(f'cd {save_dir}')

    leap_inputfile = _creat_leap_inputfile(
        prefix, ligand_file_path, prepin_file_path, frcmod_file_path, protein_file_path, save_dir)
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
