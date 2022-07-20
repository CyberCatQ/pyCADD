import os
from typing import Any

from pyCADD.utils.common import BaseFile

# If you can not import pytraj, please use the following code
# conda uninstall ambertools
# conda install -c conda-forge ambertools
import pytraj as pt
import pandas as pd
from numpy import ndarray

CPU_NUM = os.cpu_count()
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CWD = os.getcwd()


def _calc_rmsd(trajectory: Any, mask: str = '@CA', reference: Any = 0, save_dir: str = None, *args, **kwargs) -> ndarray:
    '''
    计算RMSD和RMSF

    Parameters
    ----------
    trajectory : <object pytraj.traj>
        pytraj MD轨迹对象
    mask : str
        计算RMSD的Amber mask 默认为alpha-carbon "@CA"
    reference : <object pytraj.traj> | int
        参考轨迹对象或参考轨迹索引 默认为第一帧
    save_dir : str
        保存路径 默认为当前路径

    Returns
    -------
    ndarray
        RMSD结果数组
    '''

    data_rmsd = pt.rmsd(trajectory, mask=mask, ref=reference, *args, **kwargs)

    save_dir = CWD if save_dir is None else save_dir
    output_file = os.path.join(save_dir, 'RMSD_RESULTS.csv')
    pd.DataFrame(data_rmsd, columns=['rmsd']).to_csv(
        output_file, index=False)

    return data_rmsd


def _calc_rmsf(trajectory: Any, mask: str = '@CA', options: str = 'byres', save_dir: str = None, *args, **kwargs) -> ndarray:
    '''
    计算RMSF

    Parameters
    ----------
    trajectory : <object pytraj.traj>
        pytraj MD轨迹对象
    mask : str
        计算RMSF的Amber mask 默认为alpha-carbon "@CA"
    options : str   
        计算RMSF的选项 默认为byres
    save_dir : str
        保存路径 默认为当前路径

    Returns
    -------
    ndarray
        RMSF结果数组
    '''

    data_rmsf = pt.rmsf(trajectory, mask=mask,
                        options=options, *args, **kwargs)

    save_dir = CWD if save_dir is None else save_dir
    output_file = os.path.join(save_dir, 'RMSF_RESULTS.csv')
    pd.DataFrame(data_rmsf, columns=['resnum', 'rmsf']).to_csv(
        output_file, index=False)

    return data_rmsf


def _calc_hbond(
        trajectory: Any, mask: str = ':*',
        distance: float = 3.0, angle: float = 135, options: str = None,
        save_dir: str = None, *args, **kwargs) -> tuple:
    '''
    统计轨迹氢键(H-bond)长度/角度变化

    Parameters
    ----------
    trajectory : <object pytraj.traj>
        pytraj MD轨迹对象
    mask : str
        统计轨迹氢键的Amber mask 默认为发现所有氢键
    '''
    save_dir = CWD if save_dir is None else save_dir
    avgout_file = f'{save_dir}/HBOND_RESULTS.dat'
    options = f'avgout {avgout_file} printatomnum nointramol' if options is None else options

    hbond = pt.hbond(trajectory, mask=mask, distance=distance,
                     angle=angle, options=options, *args, **kwargs)
    distance_mask, angle_mask = hbond.get_amber_mask()

    hbond_distance = _trace_distance(trajectory, distance_mask)
    hbond_angle = _trace_angle(trajectory, angle_mask)

    hbond_distance_file = os.path.join(save_dir, 'HBOND_DIS_RESULTS.csv')
    hbond_angle_file = os.path.join(save_dir, 'HBOND_ANG_RESULTS.csv')

    pd.DataFrame(hbond_distance, index=distance_mask).T.to_csv(
        hbond_distance_file, index=False)
    pd.DataFrame(hbond_angle, index=angle_mask).T.to_csv(
        hbond_angle_file, index=False)

    return hbond_distance, hbond_angle


def _trace_distance(trajectory: Any, mask: str, save: bool = False, save_dir: str = None, *args, **kwargs) -> ndarray:
    '''
    计算轨迹距离变化

    Parameters
    ----------
    trajectory : <object pytraj.traj>
        pytraj MD轨迹对象
    mask : str
        计算轨迹距离的Amber mask
    save : bool
        是否保存结果 默认为False
    save_dir : str
        保存路径 默认为当前路径

    Returns
    -------
    ndarray
        轨迹距离变化结果数组
    '''

    save_dir = CWD if save_dir is None else save_dir
    output_file = os.path.join(save_dir, f'TRACE_DIS_RESULTS_{mask}.csv')

    data_trace_distance = pt.distance(trajectory, mask=mask, *args, **kwargs)
    if save:
        pd.DataFrame(data_trace_distance, index=mask).T.to_csv(
            output_file, index=False)

    return data_trace_distance


def _trace_angle(trajectory: Any, mask: str, save: bool = False, save_dir: str = None, *args, **kwargs) -> ndarray:
    '''
    计算轨迹角度变化

    Parameters
    ----------
    trajectory : <object pytraj.traj>
        pytraj MD轨迹对象
    mask : str
        计算轨迹角度的Amber mask
    save : bool
        是否保存结果 默认为False
    save_dir : str
        保存路径 默认为当前路径

    Returns
    -------
    ndarray
        轨迹角度变化结果数组
    '''

    save_dir = CWD if save_dir is None else save_dir
    output_file = os.path.join(save_dir, 'TRACE_ANG_RESULTS.csv')

    data_trace_angle = pt.angle(trajectory, mask=mask, *args, **kwargs)
    if save:
        pd.DataFrame(data_trace_angle, index=mask).T.to_csv(
            output_file, index=False)

    return data_trace_angle


def _get_lowest_energy_info(mdout_file: BaseFile, save_dir: str = None) -> tuple:
    '''
    获取最低能量状态对应帧数

    Parameters
    ----------
    mdout_file : <object BaseFile>
        MD输出文件对象
    save_dir : str
        保存路径 默认为当前路径

    Returns
    -------
    tuple[int, int, int]
        最低能量状态对应帧数 时间 能量值
    '''

    save_dir = CWD if save_dir is None else save_dir

    os.chdir(save_dir)
    process_script_path = os.path.join(SCRIPT_DIR, 'process_mdout.perl')
    output_file = os.path.join(save_dir, 'LOWEST_ENERGY_ST_RESULTS.csv')
    os.system(f'chmod +x {process_script_path} && {process_script_path} {mdout_file.file_path}')
    
    energy_file = 'summary.EPTOT'
    energy_df = pd.read_csv(energy_file, sep='\s+', header=None, names=['time', 'energy'])
    frame = energy_df['energy'].idxmin()
    time, energy = energy_df.iloc[energy_df['energy'].idxmin(), :].values

    with open(output_file, 'w') as f:
        f.write(f'frame,time,energy\n')
        f.write(f'{frame},{time},{energy}')

    os.chdir(CWD)

    return frame, time, energy

def _extract_frame(traj: Any, frame_indices:Any, save_dir:str=None, *args, **kwargs) -> None:
    '''
    提取指定帧结构

    Parameters
    ----------
    traj : <object pytraj.traj>
        pytraj MD轨迹对象
    frame_indices : Any
        提取帧索引
    save_dir : str
        保存路径 默认为当前路径

    Returns
    -------
    None
    '''

    save_dir = CWD if save_dir is None else save_dir
    output_file = os.path.join(save_dir, 'FRAME_{frame}.pdb')
    for frame in frame_indices:
        current_output_file = output_file.format(frame=frame)
        pt.write_traj(current_output_file, traj=traj, frame_indices=[frame], overwrite=True, *args, **kwargs)
