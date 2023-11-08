import logging
import os
from datetime import datetime
from typing import Literal, Union

from pyCADD.utils.common import BaseFile
from pyCADD.utils.tool import makedirs_from_list
from pyCADD.Dynamic import core
from pyCADD.Dynamic.core import MDProcess, MinimizeProcess, NPTProcess, NVTProcess
from pyCADD.Dynamic.template import LeapInput, HeatInput, NVTInput, NPTInput, MMGBSAInput, RestrainedMinimizeInput, MinimizeInput

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

def find_amber():
    '''
    Check if AMBER is installed.
    '''
    amberhome = os.environ.get('AMBERHOME')
    _sander = os.path.exists(os.popen('which sander').read())
    _pmemd = os.path.exists(os.popen('which pmemd.cuda').read())
    _cpptraj = os.path.exists(os.popen('which cpptraj').read())
    if not all([amberhome, _sander, _pmemd, _cpptraj]):
        return False
    try:
        import pytraj as pt
    except ImportError:
        return False
    else:
        return True
    
class Processor:
    '''
    为分子动力学模拟执行预处理
    '''

    def __init__(self, apo: bool = False):
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
        if apo:
            logger.info('Initializing Processor for apo system.')
        self.apo = apo

        self.md_process_list = []

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

    def protein_prepare(self, protein_file_path: str, keep_water: bool = False) -> None:
        '''
        为动力学模拟执行蛋白质结构预处理

        Parameters
        ----------
        protein_file_path : str
            蛋白质结构文件路径
        keep_water : bool, optional
            是否保留输入结构中的水分子, 默认为 False
        '''
        protein_file = BaseFile(protein_file_path)
        self.processed_profile = core.protein_prepare(
            protein_file, save_dir=PRO_RELATED_DIR, keep_water=keep_water)
        logger.info(
            f'Protein file {self.processed_profile.file_name} has been saved in {PRO_RELATED_DIR} .')

    def molecule_prepare(
            self,
            molecule_file_path: str,
            charge: int = 0,
            multiplicity: int = 1,
            cpu_num: int = None,
            solvent: str = 'water',
            overwrite: bool = False,
            method: str = 'resp',
            keep_origin_cood: bool = False) -> None:
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
        overwrite : bool
            是否覆盖已存在的文件 默认为False
        method : str
            计算方法 默认为resp电荷
            option: resp, bcc
        keep_origin_cood : bool, optional
            是否在输出结构中保留原始坐标(计算RESP时) 
            而不使用高斯结构优化的坐标 默认为False
        '''
        cpu_num = cpu_num if cpu_num is not None else CPU_NUM
        molecule_file = BaseFile(molecule_file_path)
        if method == 'resp':
            self.processed_molfile_pdb, self.frcmod_file = core.molecule_prepare_resp2(
                molecule_file, save_dir=MOL_RELATED_DIR, charge=charge, multiplicity=multiplicity,
                cpu_num=cpu_num, solvent=solvent, overwrite=overwrite, keep_origin_cood=keep_origin_cood)
        elif method == 'bcc':
            self.processed_molfile_pdb, self.frcmod_file = core.molecule_prepare_bcc(
                molecule_file, save_dir=MOL_RELATED_DIR, charge=charge, overwrite=overwrite)
        else:
            raise ValueError(f'Method {method} is not supported.')
        logger.info(
            f'Molecule file {self.processed_molfile_pdb.file_name} has been saved in {MOL_RELATED_DIR} .')
        logger.info(
            f'Frcmod file {self.frcmod_file.file_name} has been saved in {MOL_RELATED_DIR} .')

    def load_processed_profile(self, profile_path: str) -> None:
        '''
        为Processor载入其他方法生成的蛋白结构文件

        Parameters
        ----------
        profile_path : str
            文件路径
        '''
        self.processed_profile = BaseFile(profile_path)

    def load_processed_molfile(self, molfile_path: str) -> None:
        '''
        为Processor载入其他方法生成的配体分子文件

        Parameters
        ----------
        molfile_path : str
            文件路径
        '''
        self.processed_molfile_pdb = BaseFile(molfile_path)

    def load_frcmod_file(self, file_path: str) -> None:
        '''
        为Processor载入其他方法生成的配体分子参数frcmod文件

        Parameters
        ----------
        file_path : str
            .frcmod 文件路径
        '''
        self.frcmod_file = BaseFile(file_path)

    def leap_prepare(self, prefix: str = None, box_size: float = 12.0, **kwargs) -> None:
        '''
        创建leap输入文件 并执行tleap命令

        Parameters
        ----------
        prefix : str
            leap生成文件前缀 PDBID或其他 None则为当前日期
        '''
        prefix = prefix if prefix is not None else datetime.now().strftime('%Y%m%d')
        _check_dirs = [self.processed_molfile_pdb, self.frcmod_file,
                       self.processed_profile] if not self.apo else [self.processed_profile]
        if not all(_check_dirs):
            raise RuntimeError(
                'Preparing protein or molecules before LEaP has not been done.')
        logger.info('Preparing LEaP files...')
        logger.info(f"TIP3P Water Box size: {box_size} Angstroms")
        if not self.apo:
            core.leap_prepare(
                prefix=prefix,
                ligand_file=self.processed_molfile_pdb,
                frcmod_file=self.frcmod_file,
                protein_file=self.processed_profile,
                box_size=box_size,
                save_dir=LEAP_DIR,
            )
        else:
            core.leap_prepare_for_apo(
                prefix=prefix,
                protein_file=self.processed_profile,
                box_size=box_size,
                save_dir=LEAP_DIR
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

    def get_water_resnum(self) -> list:
        '''
        获取溶剂化文件中的水分子Residue Number列表
        '''
        if self.comsolvate_pdbfile is None:
            raise ValueError(
                'Solvated complex pdb file has not been prepared/load.')
        return core._get_water_resnum(self.comsolvate_pdbfile)

    def creat_minimize_input(
            self, maxcyc: int = 10000, ncyc: int = 5000, cut: float = 8.0,
            restraint: bool = False, restraint_mask: str = None, restraint_wt: float = 2.0,
            file_name: str = None, **kwargs) -> BaseFile:
        '''
        创建能量最小化阶段输入文件

        Parameters
        ----------
        maxcyc : int, optional
            最大迭代次数, 默认10000
        ncyc : int, optional
            前ncyc步使用最速下降法, 之后使用共轭梯度法, 默认5000
        cut : float, optional
            非势能截断距离, 默认8.0
        restraint : bool, optional
            是否使用约束, 默认False
        restraint_mask : str, optional
            约束原子的Amber Mask, 启用restraint时必须指定.
        restraint_wt : float, optional
            约束力常数, 默认2.0.
        file_name : str, optional
            输入文件名, 默认minimize.in
        **kwargs : dict
            其他参数, 传入到InputFile中

        Returns
        -------
        BaseFile
            能量最小化输入文件对象
        '''
        if not restraint:
            constructor = MinimizeInput(
                maxcyc=maxcyc, ncyc=ncyc, cut=cut, **kwargs)
        else:
            if restraint_mask is None:
                raise RuntimeError(
                    'restraint_mask is required when restraint is True.')
            constructor = RestrainedMinimizeInput(
                maxcyc=maxcyc, ncyc=ncyc, cut=cut,
                restraintmask=restraint_mask, restraint_wt=restraint_wt,
                **kwargs)

        file_name = file_name if file_name is not None else 'minimize.in'
        file_path = os.path.join(INPUT_FILE_DIR, file_name)

        constructor.save(file_path)

        logger.info(
            f'Minimize process input file created for {file_name.split(".")[0]}.')
        # logger.info(
        #     f'maxcyc: {maxcyc}, ncyc: {ncyc}, restraint: {restraint}, restraint mask: {restraint_mask}')
        logger.info(f"{constructor.get_state_dict()}")

        return BaseFile(file_path)

    def creat_heat_input(
        self, tgt_temp: float = 300.0,
        heat_step: int = 9000, total_step: int = 10000,
        step_length: float = 0.002, file_name: str = None,
        restraint_wt: float = None, restraint_mask: str = None
    ) -> BaseFile:
        '''
        创建体系加热阶段输入文件

        Parameters
        ----------
        tgt_temp : float, optional
            目标温度, 默认300.0
        heat_step : int, optional
            加热步数, 默认9000
        total_step : int, optional
            总步数, 默认10000
        step_length : float, optional
            步长, 默认0.002 ps
        file_name : str, optional
            输入文件名, 默认heat.in
        restraint_wt : float, optional
            约束力常数, 默认为None, 不使用约束
        restraint_mask : str, optional
            约束原子的Amber Mask, 启用restraint时必须指定.
            
        Returns
        -------
        BaseFile
            加热阶段输入文件对象
        '''
        constructor = HeatInput(
            tgt_temperature=tgt_temp, heat_step=heat_step,
            total_step=total_step, step_length=step_length)

        file_name = file_name if file_name is not None else 'heat.in'
        file_path = os.path.join(INPUT_FILE_DIR, file_name)

        constructor.save(file_path)

        logger.info(
            f'Heat process input file created for {file_name.split(".")[0]}.')
        # logger.info(
        #     f'tgt_temp: {tgt_temp}, heat_step: {heat_step}, total_step: {total_step}')
        logger.info(f"{constructor.get_state_dict()}")

        return BaseFile(file_path)

    def creat_nvt_input(
            self, temp0: float = 300.0, total_step: int = 500000, step_length: float = 0.002,
            irest: int = 1, ntx: int = 5, file_name: str = None, **kwargs) -> BaseFile:
        '''
        创建NVT阶段输入文件

        Parameters
        ----------
        temp0 : float, optional
            初始温度, 默认300.0
        total_step : int, optional
            总步数, 默认500000, 1ns
        step_length : float, optional
            步长, 默认0.002 ps
        file_name : str, optional
            输入文件名, 默认nvt.in
        irest : int, optional
            重启标志, 默认1
        ntx : int, optional
            坐标文件输入标志, 默认5
        **kwargs : dict
            其他参数, 传入到NVT InputFile中

        Returns
        -------
        BaseFile
            NVT阶段输入文件对象
        '''
        constructor = NVTInput(
            temp0=temp0, nstlim=total_step, dt=step_length, irest=irest, ntx=ntx,
            **kwargs)

        file_name = file_name if file_name is not None else 'nvt.in'
        file_path = os.path.join(INPUT_FILE_DIR, file_name)

        constructor.save(file_path)

        logger.info(
            f'NVT process input file created for {file_name.split(".")[0]}.')
        # logger.info(
        #     f'temp0: {temp0}, total_step: {total_step}, step_length: {step_length}, irest: {irest}, ntx: {ntx}')
        logger.info(f"{constructor.get_state_dict()}")

        return BaseFile(file_path)

    def creat_npt_input(
            self, temp0: float = 300.0, total_step: int = 50000000, step_length: float = 0.002,
            irest: int = 1, ntx: int = 5, taup: float = 2.0,
            file_name: str = None, **kwargs) -> BaseFile:
        '''
        创建NPT阶段输入文件

        Parameters
        ----------
        temp0 : float, optional
            初始温度, 默认300.0
        total_step : int, optional
            总步数, 默认50000000, 100ns
        step_length : float, optional
            步长, 默认0.002 ps
        irest : int, optional
            重启标志, 默认1
        ntx : int, optional
            坐标文件输入标志, 默认5
        taup : float, optional
            压强控制时间, 默认2.0 ps
        file_name : str, optional
            输入文件名, 默认npt.in

        Returns
        -------
        BaseFile
            NPT阶段输入文件对象
        '''
        constructor = NPTInput(
            temp0=temp0, nstlim=total_step, dt=step_length, irest=irest, ntx=ntx,
            taup=taup, **kwargs)

        file_name = file_name if file_name is not None else 'npt.in'
        file_path = os.path.join(INPUT_FILE_DIR, file_name)

        constructor.save(file_path)

        logger.info(
            f'NPT process input file created for {file_name.split(".in")[0]}.')
        # logger.info(
        #     f'temp0: {temp0}, total_step: {total_step}, step_length: {step_length}, irest: {irest}, ntx: {ntx}, taup: {taup}')
        logger.info(f"{constructor.get_state_dict()}")

        return BaseFile(file_path)

    def add_process(self, input_file: Union[BaseFile, str], process_name: str = None, _type: str = None, _obj: MDProcess = None, **kwargs) -> None:
        '''
        在工作流中添加单个步骤(Minimize, NVT, NPT, etc.)

        Parameters
        ----------
        input_file : BaseFile | str
            MD步骤输入文件(或路径)
        process_name : str, optional
            工作流步骤名称, 默认为输入文件名
        _type : str, optional
            工作流步骤类型, 默认为None, 可选值为minimize, nvt, npt
        _obj : MDProcess, optional
            工作流步骤对象, 默认为None, 可选值为MinimizeProcess, NVTProcess, NPTProcess, 优先级高于_type
        **kwargs : dict
            其他参数, 写入输入文件并赋值至MDProcess属性中
        '''
        if isinstance(input_file, str):
            input_file = BaseFile(input_file)
        if _type == 'minimize':
            obj = MinimizeProcess
        elif _type == 'nvt':
            obj = NVTProcess
        elif _type == 'npt':
            obj = NPTProcess
        else:
            obj = MDProcess
        if _obj is not None:
            obj = _obj
        self.md_process_list.append(
            obj(input_file, process_name, **kwargs))
        logger.info(f'Process {process_name} has been added to Workflow.')

    def add_minimize_process(self, maxcyc: int = 10000, ncyc: int = 5000, process_name: str = 'minimize', restraint: bool = False, restraint_mask: str = None, **kwargs) -> None:
        '''
        在工作流中添加能量最小化步骤

        Parameters
        ----------
        maxcyc : int, optional
            最大迭代次数, 默认10000
        ncyc : int, optional
            前ncyc步使用最速下降法, 之后使用共轭梯度法, 默认5000
        process_name : str, optional
            工作流步骤名称, 默认为minimize.
        restraint : bool, optional
            是否添加约束, 默认为False.
        restraint_mask : str, optional
            约束原子的Amber Mask, 启用restraint时必须指定.
        **kwargs : dict
            其他参数, 传入到Minimize InputFile中
        '''
        file_name = process_name + '.in'
        self.add_process(
            self.creat_minimize_input(
                maxcyc=maxcyc, ncyc=ncyc,
                restraint=restraint, restraint_mask=restraint_mask,
                file_name=file_name, **kwargs),
            process_name,
            _type='minimize',
            **kwargs)

    def add_nvt_process(self, total_step: int = 500000, step_length: float = 0.002, process_name: str = 'nvt', is_production: bool = False, **kwargs):
        '''
        在工作流中添加NVT步骤

        Parameters
        ----------
        total_step : int, optional
            MD步骤步数, 默认50000000
        step_length : float, optional
            MD步骤步长, 默认0.002
        process_name : str, optional
            工作流步骤名称, 默认为nvt
        is_production : bool, optional
            该NVT步骤是否为生产步骤, 默认为False
        **kwargs : dict
            其他参数, 传入到NVT InputFile中
        '''
        file_name = process_name + '.in'
        self.add_process(
            self.creat_nvt_input(
                total_step=total_step, step_length=step_length, file_name=file_name, **kwargs),
            process_name,
            _type='nvt', is_production=is_production,
            **kwargs
        )

    def add_heat_process(self, tgt_temp: float = 300.0, heat_step: int = 9000, total_step: int = 10000, step_length: float = 0.002, process_name: str = 'heat', **kwargs):
        '''
        在工作流中添加加热步骤

        Parameters
        ----------
        tgt_temp : float
            加热目标温度(K), 默认300.0K
        heat_step : int, optional
            加热步数, 默认9000
        total_step : int, optional
            总加热阶段步数, 默认10000
        step_length : float, optional
            步长, 默认0.002 ps
        process_name : str, optional
            工作流步骤名称, 默认为heat
        '''
        file_name = process_name + '.in'
        self.add_process(
            self.creat_heat_input(tgt_temp=tgt_temp, heat_step=heat_step,
                                  total_step=total_step, step_length=step_length, file_name=file_name),
            process_name,
            **kwargs
        )

    def add_npt_process(self, total_step: int = 50000000, step_length: float = 0.002, process_name: str = 'npt', is_production: bool = False, **kwargs):
        '''
        在工作流中添加NPT步骤

        Parameters
        ----------
        step_num : int, optional
            步数, 默认50000000, 100ns
        step_length : float, optional
            步长, 默认0.002 ps
        process_name : str, optional
            工作流步骤名称, 默认为npt
        '''
        file_name = process_name + '.in'
        self.add_process(
            self.creat_npt_input(
                total_step=total_step, step_length=step_length, file_name=file_name, **kwargs),
            process_name, _type='npt',
            is_production=is_production,
            **kwargs
        )


class Simulator:

    def __init__(self, processor: Processor) -> None:
        self.processor = processor

        self.comsolvate_topfile = processor.comsolvate_topfile
        self.comsolvate_crdfile = processor.comsolvate_crdfile
        self.md_process_list = processor.md_process_list

        self.cuda_device = 0

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
        if len(self.md_process_list) == 0:
            raise RuntimeError('No MD process has been added.')

        if cuda_device is not None:
            self.set_cuda_device(cuda_device)
        else:
            self._apply_cuda_device()

        core._run_simulation(
            self.comsolvate_topfile, self.comsolvate_crdfile,
            self.md_process_list,
            MD_RESULT_DIR
        )

        logger.info(f'Simulation normally finished.')


class Analyzer:
    '''
    MD模拟轨迹分析器
    '''

    def __init__(
        self,
        traj_file_path: str = None, comsolvated_topfile_path: str = None, com_topfile_path: str = None,
        receptor_topfile_path: str = None, ligand_topfile_path: str = None, mdout_file_path: str = None,
    ) -> None:

        from pytraj import iterload
        self._iterload = iterload

        from pyCADD.Dynamic import analysis
        self.analyzer = analysis

        self.traj_file_path = traj_file_path
        self.traj_file = None
        self.traj = None

        self.top_file_path = comsolvated_topfile_path
        self.top_file = None
        self.top = None

        if self.top_file_path is not None and self.traj_file_path is not None:
            self.load_traj(traj_file_path, comsolvated_topfile_path)

        makedirs_from_list(self._require_dir)

        self.mdout_file_path = mdout_file_path
        self.com_topfile_path = com_topfile_path
        self.recep_topfile_path = receptor_topfile_path
        self.ligand_topfile_path = ligand_topfile_path
        
        self.apo = False if self.ligand_topfile_path is not None else True

        self.mdout_file = BaseFile(
            mdout_file_path) if mdout_file_path is not None else None
        self.com_topfile = None
        self.recep_topfile = None
        self.ligand_topfile = None

        if not self.apo:
            _check_list = [self.com_topfile_path, self.recep_topfile_path, self.ligand_topfile_path]
        else:
            _check_list = [self.com_topfile_path, self.recep_topfile_path]
        if any(_check_list):
            self.load_topfile(
                com_topfile_path=self.com_topfile_path,
                receptor_topfile_path=self.recep_topfile_path,
                ligand_topfile_path=self.ligand_topfile_path
            )

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

    def load_topfile(
        self,
        comsolvated_topfile_path: str = None, com_topfile_path: str = None,
        receptor_topfile_path: str = None, ligand_topfile_path: str = None
    ) -> None:

        if comsolvated_topfile_path is not None:
            self.top_file = BaseFile(comsolvated_topfile_path)
            logger.info(
                f'Solvated complex topology file {self.top_file.file_path} has been loaded.')
        if com_topfile_path is not None:
            self.com_topfile = BaseFile(com_topfile_path)
            logger.info(
                f'Complex topology file {self.com_topfile.file_path} has been loaded.')
        if receptor_topfile_path is not None:
            self.recep_topfile = BaseFile(receptor_topfile_path)
            logger.info(
                f'Receptor topology file {self.recep_topfile.file_path} has been loaded.')
        if ligand_topfile_path is not None:
            self.ligand_topfile = BaseFile(ligand_topfile_path)
            logger.info(
                f'Ligand topology file {self.ligand_topfile.file_path} has been loaded.')
        else:
            logger.info('Ligand topology file is not provided. Perform analysis for Apo system.')

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
        logger.info(f'Calculating RMSD...')
        logger.info(f'Amber mask: {mask}')
        logger.info(f'Reference frame: {ref}')
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

        logger.info(f'Calculating RMSF...')
        logger.info(f'Amber mask: {mask}')
        logger.info(f'Options: {options}')
        self.rmsf = self.analyzer._calc_rmsf(
            self.traj, mask=mask, options=options, save_dir=save_dir, **kwargs)

    def calc_hbond(self, mask: str = ':*', distance: float = 3.0, angle: float = 135.0, options: str = None, **kwargs) -> None:
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

        logger.info(f'Calculating and tracing H-bonds...')
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

    def extract_frame(self, frame: int, **kwargs) -> None:
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

        logger.info(f'Detecting lowest energy frame...')
        self.LE_frame, self.LE_time, self.LE_energy = self.analyzer._get_lowest_energy_info(
            mdout_file=self.mdout_file, save_dir=save_dir)

        logger.info(f'Lowest energy frame: {self.LE_frame}')
        logger.info(f'Lowest energy time: {self.LE_time}')
        logger.info(f'Lowest energy: {self.LE_energy}')
        logger.info(
            f'Extracting lowest energy structure: Frame {self.LE_frame}...')

        self.analyzer._extract_frame(self.traj, frame_indices=[
                                     self.LE_frame], save_dir=save_dir)

    def creat_energy_inputfile(
            self, start_frame: int, end_frame: int,
            job_type: Literal['free', 'entropy', 'decomp'],
            method: Literal['pb/gbsa', 'gbsa'] = None, interval: int = 10) -> None:
        '''
        创建能量计算任务输入文件

        Parameters
        ----------
        job_type : str
            能量计算任务类型
            free : 自由能计算(MM-PB/GBSA)
            entropy : 熵计算(normal mode)
            decomp : 能量分解
        method : str
            计算方法
            pb/gbsa : MM-PB/GBSA
            gbsa : MM-GBSA only
        start_frame : int
            计算分析起始帧索引号
        end_frame : int
            计算分析结束帧索引号
        interval : int
            计算分析帧间隔 默认为10
        '''
        save_dir = os.path.join(ANALYSIS_RESULT_DIR, 'energy_inputfile')
        os.makedirs(save_dir, exist_ok=True)

        if job_type == 'entropy':
            self.inputfile = core._creat_energy_inputfile(
                'nmode', startframe=start_frame,
                endframe=end_frame, interval=interval,
                save_dir=save_dir
            )
        elif job_type == 'free':
            if method is None:
                raise ValueError('Please specify method.')
            elif method == 'pb/gbsa':
                self.inputfile = core._creat_energy_inputfile(
                    'pb/gb', startframe=start_frame,
                    endframe=end_frame, interval=interval,
                    save_dir=save_dir
                )
            elif method == 'gbsa':
                self.inputfile = core._creat_energy_inputfile(
                    'gb', startframe=start_frame,
                    endframe=end_frame, interval=interval,
                    save_dir=save_dir
                )
            else:
                raise ValueError(f'Invalid method: {method}')
        elif job_type == 'decomp':
            self.inputfile = core._creat_energy_inputfile(
                'gb', startframe=start_frame,
                endframe=end_frame, interval=interval, decomp=True,
                save_dir=save_dir
            )
        else:
            raise ValueError(f'Invalid job type: {job_type}')

        return self.inputfile

    def run_energy_calc(self, output_file: str = None, decom_output_file: str = None, cpu_num: int = None, input_file: BaseFile = None) -> None:
        '''
        运行能量计算任务

        Parameters
        ----------
        output_file : str
            输出文件路径
        decom_output_file : str
            能量分解输出文件路径(仅执行能量分解有效)
        cpu_num : int
            计算使用的CPU数量
        input_file : BaseFile
            能量计算任务输入文件 默认为creat_energy_inputfile建立的文件
        '''
        save_dir = os.path.join(ANALYSIS_RESULT_DIR, 'energy_results')
        os.makedirs(save_dir, exist_ok=True)
        cpu_num = cpu_num if cpu_num is not None else CPU_NUM
        input_file = input_file if input_file is not None else self.inputfile

        if self.inputfile is None:
            raise ValueError('Please create inputfile first.')
        if self.top_file is None or self.traj_file is None:
            raise ValueError('Please load solvated complex top file first.')
        if not all([
            self.com_topfile,
            self.recep_topfile,
            self.ligand_topfile
        ]):
            raise ValueError(
                'Please load complex/receptor/ligand top file first.')

        output_file = output_file if output_file is not None else os.path.join(
            save_dir, 'FINAL_RESULTS_MMPBSA.csv')
        while True:
            i = 2
            if os.path.exists(output_file):
                output_file = os.path.join(
                    save_dir, f'{BaseFile(output_file).file_prefix}_{str(i)}.csv')
                i += 1
            else:
                break

        decom_output_file = decom_output_file if decom_output_file is not None else os.path.join(
            save_dir, 'FINAL_RESULTS_DECOMP.csv')

        core._run_energy_calculation(
            input_file=input_file, comsolvate_topfile=self.top_file,
            com_topfile=self.com_topfile, receptor_topfile=self.recep_topfile,
            ligand_topfile=self.ligand_topfile, traj_file=self.traj_file,
            output_filepath=output_file, decom_output_filepath=decom_output_file,
            cpu_num=cpu_num
        )
        logger.info(f'Calculating normally finished.')
