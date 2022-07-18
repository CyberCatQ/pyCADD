import logging
import os
import time
from datetime import datetime

from pyCADD.Dynamic import core
from pyCADD.utils.common import BaseFile
from pyCADD.utils.tool import makedirs_from_list

logger = logging.getLogger('pyCADD.Dynamic')

CWD = os.getcwd()
CPU_NUM = os.cpu_count()
SANDER = 'pmemd.cuda'

PRO_RELATED_DIR = os.path.join(CWD, 'protein')
MOL_RELATED_DIR = os.path.join(CWD, 'molecule')
LEAP_DIR = os.path.join(CWD, 'leap')
INPUT_FILE_DIR = os.path.join(CWD, 'input_file')
MD_RESULT_DIR = os.path.join(CWD, 'md_result')

ANALYSIS_RESULT_DIR = os.path.join(CWD, 'md_analysis')


class Processor:
    '''
    为分子动力学模拟执行预处理
    '''

    def __init__(self):
        '''
        初始化
        '''
        makedirs_from_list(self._required_dirs)

        self.processed_profile = None
        self.processed_molfile_pdb = None
        self.frcmod_file = None
        self.comsolvate_crdfile = None
        self.comsolvate_topfile = None
        self.comsolvate_pdbfile = None

    @property
    def _required_dirs(self):
        return [
            PRO_RELATED_DIR,
            MOL_RELATED_DIR,
            LEAP_DIR,
            INPUT_FILE_DIR,
            MD_RESULT_DIR
        ]

    @staticmethod
    def cover_to_pdb(file_path: str) -> None:
        '''
        将mae结构文件转换为pdb格式

        Parameters
        ----------
        file_path : str
            mae结构文件路径
        '''
        core._convert_mae_to_pdb(file_path)

    def protein_prepare(self, pretein_file_path: str) -> None:
        '''
        为动力学模拟执行蛋白质结构预处理

        Parameters
        ----------
        pretein_file_path : str
            蛋白质结构文件路径
        '''
        protein_file = BaseFile(pretein_file_path)
        self.processed_profile = core.protein_prepare(
            protein_file, save_dir=PRO_RELATED_DIR)
        logger.info(
            f'Protein file {self.processed_profile.file_name} has been saved in {PRO_RELATED_DIR} .')

    def molecule_prepare(self, molecule_file_path: str, charge: int = 0, multiplicity: int = 1, cpu_num: int = None, solvent: str = 'water') -> None:
        '''
        为动力学模拟执行小分子结构预处理

        Parameters
        ----------
        molecule_file_path : str
            分子结构文件路径
        charge : int
            电荷数 默认为0
        multiplicity : int
            自旋多重度 默认为1
        cpu_num : int
            计算核数 默认为CPU核数
        solvent : str
            计算RESP电荷时液相的溶剂 默认为water
        '''
        cpu_num = cpu_num if cpu_num is not None else CPU_NUM
        molecule_file = BaseFile(molecule_file_path)
        self.processed_molfile_pdb, self.frcmod_file = core.molecule_prepare(
            molecule_file, save_dir=MOL_RELATED_DIR, charge=charge, multiplicity=multiplicity,
            cpu_num=cpu_num, solvent=solvent)
        logger.info(
            f'Molecule file {self.processed_molfile_pdb.file_name} has been saved in {MOL_RELATED_DIR} .')
        logger.info(
            f'Frcmod file {self.frcmod_file.file_name} has been saved in {MOL_RELATED_DIR} .')

    def leap_prepare(self, prefix: str = None) -> None:
        '''
        创建leap输入文件 并执行tleap命令

        Parameters
        ----------
        prefix : str
            leap生成文件前缀 PDBID或其他 None则为当前日期
        '''
        prefix = prefix if prefix is not None else datetime.now().strftime('%Y%m%d')
        if not all([
            self.processed_molfile_pdb,
            self.frcmod_file,
            self.processed_profile
        ]):
            raise RuntimeError(
                'Preparing protein or molecules before LEaP has not been done.')

        core.leap_prepare(
            prefix,
            self.processed_molfile_pdb,
            self.frcmod_file,
            self.processed_profile,
            LEAP_DIR
        )

        self.comsolvate_pdbfile = BaseFile(
            os.path.join(LEAP_DIR, f'{prefix}_comsolvate.pdb'))
        self.comsolvate_topfile = BaseFile(
            os.path.join(LEAP_DIR, f'{prefix}_comsolvate.prmtop'))
        self.comsolvate_crdfile = BaseFile(
            os.path.join(LEAP_DIR, f'{prefix}_comsolvate.inpcrd'))

        logger.info(f'LEaP files have been saved in {LEAP_DIR} .')

    def _set_prepared_file(self, file_path: str, file_type: str) -> None:
        '''
        直接设定已准备好的文件

        Parameters
        ----------
        file_path : str
            文件路径
        file_type : str
            文件类型
            protein: 蛋白质PDB结构文件(_leap.pdb)
            molecule: 小分子PDB结构文件(_out.pdb)
            frcmod: amber参数文件(.frcmod)
            comsolvate_pdb: 含水蛋白结构文件(_comsolvate.pdb)
            comsolvate_top: 含水蛋白拓扑文件(_comsolvate.prmtop)
            comsolvate_crd: 含水蛋白坐标文件(_comsolvate.inpcrd)
        '''
        _file = BaseFile(file_path)
        if file_type == 'protein':
            self.processed_profile = _file
        elif file_type == 'molecule':
            self.processed_molfile_pdb = _file
        elif file_type == 'frcmod':
            self.frcmod_file = _file
        elif file_type == 'comsolvate_pdb':
            self.comsolvate_pdbfile = _file
        elif file_type == 'comsolvate_top':
            self.comsolvate_topfile = _file
        elif file_type == 'comsolvate_crd':
            self.comsolvate_crdfile = _file
        else:
            raise RuntimeError(f'{file_type} is not a valid file type.')

    def set_comsolvate_file(self, file_path: str, file_type: str) -> None:
        '''
        设定含水蛋白结构相关文件

        Parameters
        ----------
        file_path : str
            文件路径
        file_type : str
            文件类型
            pdb: 含水蛋白结构文件(_comsolvate.pdb)
            top: 含水蛋白拓扑文件(_comsolvate.prmtop)
            crd: 含水蛋白坐标文件(_comsolvate.inpcrd)
        '''
        if file_type == 'pdb':
            self._set_prepared_file(file_path, 'comsolvate_pdb')
            logger.info(f'Set comsolvate pdb file: {file_path}')
        elif file_type == 'top':
            self._set_prepared_file(file_path, 'comsolvate_top')
            logger.info(f'Set comsolvate top file: {file_path}')
        elif file_type == 'crd':
            self._set_prepared_file(file_path, 'comsolvate_crd')
            logger.info(f'Set comsolvate crd file: {file_path}')
        else:
            self._set_prepared_file(file_path, file_type)
            logger.info(f'Set {file_type} file: {file_path}')


class Simulator:

    def __init__(self, processor: Processor) -> None:
        self.processor = processor

        self.comsolvate_pdbfile = processor.comsolvate_pdbfile
        self.comsolvate_topfile = processor.comsolvate_topfile
        self.comsolvate_crdfile = processor.comsolvate_crdfile

        self.step_a_inputfile = None
        self.step_b_inputfile = None
        self.step_c_inputfile = None
        self.step_nvt_inputfile = None
        self.step_npt_inputfile = None

        self.cuda_device = 0

    def creat_input_file(self, step_num: int = 50000000, step_length: float = 0.002):
        '''
        创建动力学模拟所需的各类输入文件

        Parameters
        ----------
        step_num : int
            模拟步数
        step_length : float
            模拟步长

        Notes
        -----
        模拟总时长 = step_num * step_length ps
        默认100ns
        '''
        water_resnum = core._get_water_resnum(self.comsolvate_pdbfile)
        (self.step_a_inputfile, self.setp_b_inputfile,
         self.setp_c_inputfile, self.step_nvt_inputfile,
         self.step_npt_inputfile) = core._creat_md_inputfile(
            water_resnum, step_num,
            step_length, INPUT_FILE_DIR)

        logger.info(f'Input files have been saved in {INPUT_FILE_DIR}.')

    def shwo_cuda_device(self) -> None:
        '''
        显示GPU信息
        '''

        gpu_info = os.popen('nvidia-smi').read()
        print(gpu_info)

    def _apply_cuda_device(self) -> None:
        '''
        设定模拟使用的GPU设备

        Parameters
        ----------
        cuda_device : int
            GPU设备编号
        '''
        cuda_device = self.cuda_device
        print(os.popen(f'nvidia-smi -i {cuda_device}').read())
        os.environ['CUDA_VISIBLE_DEVICES'] = str(cuda_device)
        logger.info(f'Using GPU device {cuda_device}')

    def set_cuda_device(self, cuda_device: int) -> None:
        '''
        设定模拟使用的GPU设备

        Parameters
        ----------
        cuda_device : int
            GPU设备编号
        '''
        self.cuda_device = cuda_device
        logger.info(f'Set GPU device: {cuda_device}')
        self._apply_cuda_device()

    def run_simulation(self, cuda_device: int = None) -> None:
        '''
        启动分子动力学模拟流程

        Parameters
        ----------
        cuda_device : int
            GPU设备编号 默认为0
        '''
        if not all([
            self.step_a_inputfile,
            self.setp_b_inputfile,
            self.setp_c_inputfile,
            self.step_nvt_inputfile,
            self.step_npt_inputfile
        ]):
            raise RuntimeError('Creating input files has not been done.')

        if cuda_device is not None:
            self.set_cuda_device(cuda_device)
        else:
            self._apply_cuda_device()

        start_time = time.time()
        core._run_simulation(
            self.comsolvate_topfile, self.comsolvate_crdfile,
            self.step_a_inputfile, self.setp_b_inputfile,
            self.setp_c_inputfile, self.step_nvt_inputfile,
            self.step_npt_inputfile,
            MD_RESULT_DIR
        )
        end_time = time.time()

        duration = time.localtime(end_time - start_time)

        logger.info(
            f'Start: {time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(start_time))}')
        logger.info(
            f'End: {time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(end_time))}')
        logger.info(f'Simulation time: {time.strftime("%H:%M:%S", duration)}')
        logger.info(f'Simulation normally finished.')


class Analyzer:
    '''
    MD模拟轨迹分析器
    '''

    def __init__(self) -> None:

        from pytraj import iterload
        self._iterload = iterload

        from pyCADD.Dynamic import analysis
        self.analyzer = analysis

        self.traj_file_path = None
        self.traj_file = None
        self.traj = None

        self.top_file_path = None
        self.top_file = None
        self.top = None
        makedirs_from_list(self._require_dir)

        self.mdout_file = None

        self.rmsd = None
        self.rmsf = None
        self.hbond = None
        self.distance = None
        self.angle = None

    @property
    def _require_dir(self):
        return [
            ANALYSIS_RESULT_DIR
        ]

    def load_traj(self, traj_file_path: str, top_file_path: str) -> None:
        '''
        加载轨迹文件

        Parameters
        ----------
        traj_file_path : str
            轨迹文件路径
        top_file_path : str
            拓扑文件路径
        '''

        self.traj_file = BaseFile(traj_file_path)
        self.top_file = BaseFile(top_file_path)

        self.traj = self._iterload(
            self.traj_file.file_path, self.top_file.file_path)

        logger.info(
            f'Trajectory file {self.traj_file.file_path} has been loaded.')
        logger.info(
            f'Topology file {self.top_file.file_path} has been loaded.')
        print('Trajectory Info:\n', self.traj)
    
    def load_mdout(self, mdout_file_path: str) -> None:
        '''
        加载MD模拟输出文件

        Parameters
        ----------
        mdout_file_path : str
            MD模拟输出文件路径
        '''

        self.mdout_file = BaseFile(mdout_file_path)
        logger.info(
            f'MD output file {self.mdout_file.file_path} has been loaded.')



    def calc_rmsd(self, mask: str = '@CA', ref: int = 0, **kwargs) -> None:
        '''
        计算RMSD

        Parameters
        ----------
        mask : str
            计算RMSD的Amber mask
        ref : int
            参考轨迹索引号 默认第一帧0
        '''

        save_dir = os.path.join(ANALYSIS_RESULT_DIR, 'rmsd')
        os.makedirs(save_dir, exist_ok=True)
        self.rmsd = self.analyzer._calc_rmsd(
            self.traj, mask=mask, reference=ref, save_dir=save_dir, **kwargs)

    def calc_rmsf(self, mask: str = '@CA', options: str = 'byres', **kwargs) -> None:
        '''
        计算RMSF

        Parameters
        ----------
        mask : str
            计算RMSF的Amber mask
        options : str
            RMSF计算选项 默认byres
        '''
        save_dir = os.path.join(ANALYSIS_RESULT_DIR, 'rmsf')
        os.makedirs(save_dir, exist_ok=True)
        self.rmsf = self.analyzer._calc_rmsf(
            self.traj, mask=mask, options=options, save_dir=save_dir, **kwargs)

    def calc_hbond(self, mask: str = None, distance: float = 3.0, angle: float = 135.0, options: str = None, **kwargs) -> None:
        '''
        检测并追踪 Hbond 键长、键角

        Parameters
        ----------
        mask : str
            指定原子间的Amber mask 默认为自动识别 cutoff < distance 的原子
        options : str
            Hbond计算选项 默认为
                avgout OUTPUTFILE 输出氢键平均信息文件   
                printatomnum 打印原子序号   
                nointramol 仅计算分子间氢键   
        '''
        save_dir = os.path.join(ANALYSIS_RESULT_DIR, 'hbond')
        os.makedirs(save_dir, exist_ok=True)
        self.hbond = self.analyzer._calc_hbond(
            self.traj, mask=mask, distance=distance, angle=angle, options=options, save_dir=save_dir, **kwargs)

    def trace_distance(self, mask: str, **kwargs) -> None:
        '''
        追踪指定原子之间的距离

        Parameters
        ----------
        mask : str
            指定原子Amber mask
        '''
        save_dir = os.path.join(ANALYSIS_RESULT_DIR, 'distance')
        os.makedirs(save_dir, exist_ok=True)
        self.distance = self.analyzer._trace_distance(
            self.traj, mask=mask, save=True, save_dir=save_dir, **kwargs)

    def trace_angle(self, mask: str, **kwargs) -> None:
        '''
        追踪指定原子之间的键角

        Parameters
        ----------
        mask : str
            指定原子Amber mask
        '''
        save_dir = os.path.join(ANALYSIS_RESULT_DIR, 'angle')
        os.makedirs(save_dir, exist_ok=True)
        self.angle = self.analyzer._trace_angle(
            self.traj, mask=mask, save=True, save_dir=save_dir, **kwargs)

    def extract_frame(self, frame:int, **kwargs) -> None:
        '''
        提取指定帧
        帧索引开始于0 结束于最大帧数量-1

        Parameters
        ----------
        frame : int
            提取帧索引号
        '''
        save_dir = os.path.join(ANALYSIS_RESULT_DIR, 'frame_structures')
        os.makedirs(save_dir, exist_ok=True)
        self.frame = self.analyzer._extract_frame(
            traj=self.traj, frame_indices=[frame], save_dir=save_dir, **kwargs)

    def extract_frames(self, start: int, end: int, **kwargs) -> None:
        '''
        提取指定帧范围
        帧索引开始于0 结束于最大帧数量-1

        Parameters
        ----------
        start : int
            开始帧索引号
        end : int
            结束帧索引号
        '''
        save_dir = os.path.join(ANALYSIS_RESULT_DIR, 'frame_structures')
        os.makedirs(save_dir, exist_ok=True)
        indices = range(start, end + 1)
        self.frames = self.analyzer._extract_frame(
            traj=self.traj, frame_indices=indices, save_dir=save_dir, **kwargs)
    
    def extract_lowest_energy_st(self) -> None:
        '''
        提取最低能量结构
        '''
        save_dir = os.path.join(ANALYSIS_RESULT_DIR, 'LE_structures')
        os.makedirs(save_dir, exist_ok=True)
        if self.mdout_file is None:
            raise ValueError('Please load mdout file first.')

        self.LE_frame, self.LE_time, self.LE_energy = self.analyzer._get_lowest_energy_info(
            mdout_file=self.mdout_file, save_dir=save_dir)
        
        self.analyzer._extract_frame(self.traj, frame_indices=[self.LE_frame], save_dir=save_dir)