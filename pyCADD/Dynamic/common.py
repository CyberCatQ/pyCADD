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
        self.prepin_file = None
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
        logger.info(f'Protein file {self.processed_profile.file_name} has been saved in {PRO_RELATED_DIR}.')

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
        '''
        cpu_num = cpu_num if cpu_num is not None else CPU_NUM
        molecule_file = BaseFile(molecule_file_path)
        self.processed_molfile_pdb, self.prepin_file, self.frcmod_file = core.molecule_prepare(
            molecule_file, save_dir=MOL_RELATED_DIR, charge=charge, multiplicity=multiplicity,
            cpu_num=cpu_num, solvent=solvent)
        logger.info(f'Molecule file {self.processed_molfile_pdb.file_name} has been saved in {MOL_RELATED_DIR}.')
        logger.info(f'Prepin file {self.prepin_file.file_name} has been saved in {MOL_RELATED_DIR}.')
        logger.info(f'Frcmod file {self.frcmod_file.file_name} has been saved in {MOL_RELATED_DIR}.')

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
            self.prepin_file,
            self.frcmod_file,
            self.processed_profile
        ]):
            raise RuntimeError(
                'Preparing protein or molecules before LEaP has not been done.')

        core.leap_prepare(
            prefix,
            self.processed_molfile_pdb,
            self.prepin_file,
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
        
        logger.info(f'LEaP files have been saved in {LEAP_DIR}.')

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
            prepin: amber参数文件(.prepin)
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
        elif file_type == 'prepin':
            self.prepin_file = _file
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

        logger.info(f'Start: {time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(start_time))}')
        logger.info(f'End: {time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(end_time))}')
        logger.info(f'Simulation time: {time.strftime("%H:%M:%S", duration)}')
        logger.info(f'Simulation normally finished.')
