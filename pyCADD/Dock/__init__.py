import logging
import os

from pyCADD.Dock.common import LigandFile, MultiInputFile, PDBFile, get_input_pdbid
from pyCADD.Dock.ensemble import _Console
from pyCADD.Dock.core import minimize, grid_generate, dock, calc_mmgbsa, calc_admet, calc_volume
from pyCADD.Dock.data import extra_docking_data, extra_admet_data, save_docking_data, save_admet_data
from pyCADD.utils.tool import download_pdb

logger = logging.getLogger('pyCADD.Dock.Docker')

class Docker:
    '''
    Ligand Docking 控制台对象
    '''
    
    def __init__(self, pdbid:str=None) -> None:
        self._pdbid = pdbid
        self._pdb_file = None
        self._lig_info = None

        self.precision = None                         # 默认对接精度SP
        self.calc_rmsd = False                        # 默认不计算RMSD
        self.minimized_file = None                    # Minimized化完成的文件名
        self.grid_file = None                         # 格点文件名
        self.lig_file = None                          # 内源配体文件名
        self.recep_file = None                        # 受体文件名
        self.dock_file = None                         # 对接结果文件名
        self.mmgbsa_file = None                       # 结合能计算结果文件名
        self.sitemap_file = None                      # 结合口袋体积计算结果文件名
        self.admet_file = None                        # ADMET计算结果文件名
        self.data_dic = None                          # 一般计算结果字典
        self.admet_dic = None                         # ADMET计算结果字典

    @property
    def pdbid(self):
        if self._pdbid is None:
            self._pdbid = get_input_pdbid()
        return self._pdbid
    
    @property
    def pdb_file_path(self):
        return os.path.join(os.getcwd(), self.pdbid + '.pdb')

    def download_pdb(self) -> None:
        '''
        下载PDB文件
        '''
        logger.info(f'Prepare to download pdb: {self.pdbid}')
        if os.path.exists(self.pdb_file_path):
            logger.info(f'{self.pdbid} PDB file exists.')
            self._pdb_file = PDBFile(self.pdb_file_path)
            return self._pdb_file

        download_pdb(self.pdbid, os.getcwd())
        self._pdb_file = PDBFile(self.pdb_file_path)
        logger.info(f'Download {self.pdbid} pdb done.')
        return self._pdb_file

    def _get_lig_info(self) -> None:
        '''
        获取配体信息
        '''
        pdbid = self.pdbid
        ligand_id = input('Enter ligand ID: ')
        self._lig_info = self.pdb_file.get_lig(ligand_id)
        logger.info(f'Current Ligand: {self.pdb_file.ligid} {self.pdb_file.lig_resnum}')
        return self._lig_info

    @property
    def pdb_file(self):
        return self._pdb_file if self._pdb_file is not None else self.download_pdb()

    @property
    def lig_info(self):
        return self._lig_info if self._lig_info is not None else self._get_lig_info()

    @property
    def lig_name(self):
        return self.lig_info['id']

    def set_precision(self, precision: str) -> None:
        '''
        设置对接精度
        '''
        logger.info(f'Set precision: {precision}')
        self.precision = precision
    
    def set_calc_rmsd(self, calc_rmsd: bool) -> None:
        '''
        设置是否计算RMSD
        '''
        logger.info(f'Set calc_rmsd: {calc_rmsd}')
        self.calc_rmsd = calc_rmsd

    def keep_chain(self, *args, **kwargs) -> None:
        '''
        删改PDB晶体结构 仅保留配体所在链

        Parameters
        ----------
        *args : list, optional
            链名, 默认None
        **kwargs : dict, optional
            链名, 默认None
        '''
        _chain = self.lig_info['chain']
        logger.info(f'Prepare to keep crystal {self.pdbid} chain {_chain}')
        self._pdb_file = self.pdb_file.keep_chain(_chain)
        logger.info(f'{self.pdbid} Keep chain done.')

    def minimize(self, side_chain:bool=True, missing_loop:bool=True, del_water:bool=True, *args, **kwargs) -> None:
        '''
        优化晶体并执行能量最小化

        Parameters
        ----------
        side_chain : bool, optional
            是否优化侧链, 默认True
        missing_loop : bool, optional
            是否优化缺失的loop, 默认True
        del_water : bool, optional
            是否删除水分子, 默认True
        *args : list, optional
            优化参数, 默认None
        **kwargs : dict, optional
            优化参数, 默认None
        '''
        logger.info(f'Prepare to optimize structure: {self.pdbid}')
        self.minimized_file = minimize(self.pdb_file, side_chain=side_chain, missing_loop=missing_loop, del_water=del_water, *args, **kwargs)
        logger.info(f'{self.pdbid} Minimized structure done.')
    
    def grid_generate(self, *args, **kwargs) -> None:
        '''
        生成格点

        Parameters
        ----------
        *args : list, optional
            格点参数, 默认None
        **kwargs : dict, optional
            格点参数, 默认None
        '''
        logger.info(f'Prepare to generate grid: {self.pdbid}')
        self.minimized_file.ligid = self.lig_name
        self.minimized_file.lig_resnum = self.lig_info['resid']
        self.grid_file = grid_generate(self.minimized_file, *args, **kwargs)
        logger.info(f'{self.pdbid} Grid generated.')
    
    def split_complex(self) -> None:
        '''
        拆分复合物结构
        '''
        logger.info(f'Prepare to split complex: {self.pdbid}')
        self.recep_file, self.lig_file = self.minimized_file.split(self.lig_name)
        logger.info(f'{self.pdbid} Split complex done.')
    
    def dock(self, docking_ligand:LigandFile=None, *args, **kwargs) -> None:
        '''
        对接

        Parameters
        ----------
        docking_ligand : str, optional
            对接配体, 默认None将执行共结晶配体的回顾性对接
        *args : list, optional
            对接参数, 默认None
        **kwargs : dict, optional
            对接参数, 默认None
        '''
        docking_ligand = docking_ligand if docking_ligand is not None else self.lig_file
        logger.info(f'Prepare to dock: {self.pdbid}')
        logger.info(f'Docking ligand: {docking_ligand.ligand_name}')
        self.dock_file = dock(self.grid_file, self.lig_file, self.precision, self.calc_rmsd, *args, **kwargs)
        logger.info(f'{self.pdbid}-{docking_ligand.ligand_name}-{self.precision} Docking done.')
    
    def calc_mmgbsa(self, *args, **kwargs) -> None:
        '''
        计算结合能

        Parameters
        ----------
        *args : list, optional
            计算结合能参数, 默认None
        **kwargs : dict, optional
            计算结合能参数, 默认None
        '''
        logger.info(f'Prepare to calculate MMGBSA: {self.pdbid}')
        self.mmgbsa_file = calc_mmgbsa(self.dock_file, *args, **kwargs)
        logger.info(f'{self.pdbid} MMGBSA calculated.')
    
    def calc_volume(self, *args, **kwargs) -> None:
        '''
        计算结合口袋体积

        Parameters
        ----------
        *args : list, optional
            计算结合口袋体积参数, 默认None
        **kwargs : dict, optional
            计算结合口袋体积参数, 默认None
        '''
        logger.info(f'Prepare to calculate volume: {self.pdbid}')
        self.sitemap_file = calc_volume(self.recep_file, self.lig_file, *args, **kwargs)
        logger.info(f'{self.pdbid} Volume calculated.')
    
    def calc_admet(self, ligand_file:LigandFile=None, *args, **kwargs) -> None:
        '''
        计算ADMET

        Parameters
        ----------
        ligand_file : LigandFile
            计算ADMET预测结果的配体文件 默认为当前共结晶配体
        *args : list, optional
            计算ADMET参数, 默认None
        **kwargs : dict, optional
            计算ADMET参数, 默认None
        '''
        ligand_file = ligand_file if ligand_file is not None else self.lig_file
        logger.info(f'Prepare to calculate ADMET: {self.lig_file.file_name}')
        self.admet_file = calc_admet(ligand_file, *args, **kwargs)
        logger.info(f'{self.lig_file.file_name} ADMET calculated.')

    def extra_docking_data(self, *args, **kwargs) -> None:
        '''
        提取对接数据

        Parameters
        ----------
        *args : list, optional
            提取对接数据参数, 默认None
        **kwargs : dict, optional
            提取对接数据参数, 默认None
        '''
        logger.info(f'Prepare to extract docking data: {self.dock_file.file_name}')
        self.docking_data = extra_docking_data(self.dock_file, *args, **kwargs)
        logger.info(f'{self.dock_file.file_name} Docking data extracted.')
        return self.docking_data
    
    def extra_admet_data(self, *args, **kwargs) -> None:
        '''
        提取ADMET数据

        Parameters
        ----------
        *args : list, optional
            提取ADMET数据参数, 默认None
        **kwargs : dict, optional
            提取ADMET数据参数, 默认None
        '''
        logger.info(f'Prepare to extract ADMET data: {self.admet_file.file_name}')
        self.admet_data = extra_admet_data(self.admet_file, *args, **kwargs)
        logger.info(f'{self.admet_file.file_name} ADMET data extracted.')
        return self.admet_data

    def save_docking_data(self, *args, **kwargs) -> None:
        '''
        保存对接数据

        Parameters
        ----------
        *args : list, optional
            保存对接数据参数, 默认None
        **kwargs : dict, optional
            保存对接数据参数, 默认None
        '''
        logger.info(f'Prepare to save docking data: {self.dock_file.file_name}')
        output = save_docking_data(self.dock_file, *args, **kwargs)
        logger.info(f'{output} Docking data saved.')
    
    def save_admet_data(self, *args, **kwargs) -> None:
        '''
        保存ADMET数据

        Parameters
        ----------
        *args : list, optional
            保存ADMET数据参数, 默认None
        **kwargs : dict, optional
            保存ADMET数据参数, 默认None
        '''
        logger.info(f'Prepare to save ADMET data: {self.admet_file.file_name}')
        output = save_admet_data(self.admet_file, *args, **kwargs)
        logger.info(f'{output} ADMET data saved.')


class MultiDocker(_Console):
    def __init__(self, input_file:MultiInputFile, *args, **kwargs):
        super().__init__(input_file, *args, **kwargs)
