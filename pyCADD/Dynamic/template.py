# This file is a database of all the templates used in the Dynamic module.
# using str.format() to replace the placeholders with the actual values.

from abc import ABC, abstractmethod
from typing import Literal, Union

#########################
# Basic Class Definition
#########################


class BaseConstructor(ABC):
    def __init__(self, state_dict: dict) -> None:
        """
        Base class for all the constructors.
        """
        self._state_dict = state_dict.copy()
        self._state_keys = list(self._state_dict.keys())
        for key in self._state_keys:
            self.__dict__[key] = self._state_dict[key]

    def get_state_dict(self) -> dict:
        """
        Return the current state of the constructor as a dict.

        Returns
        -------
        dict
            Current state of the constructor.
        """
        return {key: self.__dict__[key] for key in self._state_keys}

    @abstractmethod
    def to_string(self) -> str:
        """
        Return the current state with formated template as a string.

        Returns
        -------
        str
            Formated template as a string.
        """
        ...

    def to_str(self) -> str:
        """
        Same as to_string().
        """
        return self.to_string()

    def to_dict(self) -> dict:
        """
        Return the current state as a dict with current state.

        Returns
        -------
        dict
            Current state as a dict.
        """
        return self.get_state_dict()

    def save(self, file_path: str) -> None:
        """
        Save the current state as a file.
        """
        with open(file_path, "w") as f:
            f.write(self.to_string())

    def add_attr(self, **kwargs) -> None:
        """
        Add attributes to the constructor.
        """
        for key in kwargs.keys():
            if key in self._state_keys:
                raise KeyError(
                    f"Attribute {key} already exists in the current state.")
        self.__dict__.update(kwargs)
        self._state_keys.extend(list(kwargs.keys()))

    def get_attr(self, attr_name: str) -> str:
        """
        Get the value of the attribute.
        """
        return self.get_state_dict()[attr_name]

    def set_attr(self, **kwargs) -> None:
        """
        Set the value of the attribute.
        """
        for key in kwargs.keys:
            if key not in self._state_keys:
                raise KeyError(
                    f"Attribute {key} can not be found in the current state.")
        self.__dict__.update(kwargs)

    def del_attr(self, attr_name: str) -> None:
        """
        Delete the attribute from the current state.
        """
        if attr_name not in self._state_keys:
            raise KeyError(
                f"Attribute {attr_name} can not be found in the current state.")
        del self.__dict__[attr_name]
        self._state_keys.remove(attr_name)

_LEAP_BASE_TEMPLATE_APO = """

############Protein###########
pro = loadpdb {protein_file_path}
saveamberparm pro {file_prefix}_pro.prmtop {file_prefix}_pro.inpcrd
##############################

##########Complex#############
com = combine {pro_lig}
savepdb com {file_prefix}_com.pdb
saveamberparm com {file_prefix}_com.prmtop {file_prefix}_com.inpcrd
solvatebox com {box_type} {box_size}
check com
addions com {add_ions_type} 0
saveamberparm com {file_prefix}_comsolvate.prmtop {file_prefix}_comsolvate.inpcrd
##############################

quit

"""
_LEAP_BASE_TEMPLATE_LIGAND = """
############Ligand############
loadamberparams {frcmod_file_path}
lig = loadmol2 {ligand_file_path}
saveamberparm lig {file_prefix}_lig.prmtop {file_prefix}_lig.inpcrd
##############################
"""
_LEAP_BASE_TEMPLATE = _LEAP_BASE_TEMPLATE_LIGAND + _LEAP_BASE_TEMPLATE_APO


class LeapConstructor(BaseConstructor):
    def __init__(self, state_dict: dict) -> None:
        """
        LEaP 文件用于生成蛋白-配体复合物的拓扑和坐标文件。
        LEaP file is used to generate the topology and coordinate files for the protein-ligand complex.
        state_dict = {
            "force_field": [...] ,
            "frcmod_file_path": ... ,
            "ligand_file_path": ... ,
            "file_prefix": ... ,
            "protein_file_path": ... ,
            "box_type": ... ,
            "box_size": ... ,
            "add_ions_type": ... 
            }
        """
        super().__init__(state_dict)
        self.add_attr(pro_lig=r"{pro lig}")

    @property
    def base_template(self) -> str:
        """
        Return the base template of the LEaP file.
        """
        if self.get_state_dict().get("ligand_file_path") is None:
            return _LEAP_BASE_TEMPLATE_APO
        else:
            return _LEAP_BASE_TEMPLATE

    def _get_force_field_template(self) -> str:
        """
        Return the force field template.
        """
        _tmplate = "\n".join(
            [f"source leaprc.{ff}" for ff in self.get_state_dict()["force_field"]])
        return _tmplate

    def to_string(self) -> str:
        output_str = self._get_force_field_template()
        output_str += self.base_template.format(**self.get_state_dict())
        return output_str


class SimulationConstructor(BaseConstructor):

    def __init__(self, state_dict: dict, type: Literal['cntrl', 'wt'], end_cfg: bool = True, title: str = None) -> None:
        """
        MD 模拟输入文件用于执行 MD 模拟所需的能量最小化、升温加热、平衡及正式模拟。
        MD Simulation input file is used to perform the energy minimization, heating, equilibration and production.

        Parameters
        ----------
        state_dict : dict
            The state dict of the simulation stage.
        type : Literal['cntrl', 'wt']
            The type of the simulation stage.
        end_cfg: bool, optional
            Whether to add string "END" to the end of input file. Defaults to True.
        title: str, optional
            The title of the simulation stage. Defaults to None.
        """
        super().__init__(state_dict)
        self._type = type
        self._end_cfg = end_cfg
        self.title = title

    def to_string(self) -> str:
        if self.title is not None:
            output_str = f"{self.title}\n&{self._type}\n"
        else:
            output_str = f"&{self._type}\n"
        output_str += ",\n".join([f"{str(key)}={str(value)}" for key,
                                 value in self.get_state_dict().items()])
        output_str += "\n/\n"
        if self._end_cfg:
            output_str += '\nEND'
        return output_str


class AnalysisConstructor(BaseConstructor):
    def __init__(self, state_dict: dict, type: Literal['general', 'gb', 'pb', 'decomp', 'nmode'], title: str = None) -> None:
        """
        MD 轨迹分析输入文件用于分析 MD 轨迹。
        MD Trajectory Analysis input file is used to analyze the MD trajectory.
        """
        super().__init__(state_dict)
        self._type = type
        self.title = title

    def to_string(self) -> str:
        output_str = f"# {self.title}\n" if self.title is not None else ""
        output_str += f"&{self._type}\n"
        output_str += ",\n".join([f" {str(key)}={str(value)}" for key,
                                 value in self.get_state_dict().items()])
        output_str += '\n/\n'
        return output_str


class MultiConstructorManager:
    def __init__(self, constructor_list: list = None) -> None:
        self.constructor_list = constructor_list if constructor_list is not None else []

    @property
    def _state_dict(self) -> dict:
        _state_dict = dict()
        if len(self.constructor_list) == 0:
            return _state_dict
        else:
            for constructor in self.constructor_list:
                _state_dict.update(constructor.get_state_dict())
            return _state_dict
    
    def is_empty(self) -> bool:
        return len(self.constructor_list) == 0

    def add_constructor(self, constructor: BaseConstructor) -> None:
        self.constructor_list.append(constructor)

    def del_constructor(self, constructor: BaseConstructor) -> None:
        self.constructor_list.remove(constructor)

    def to_string(self) -> str:
        output_str = "".join([tmp.to_string()
                             for tmp in self.constructor_list])
        return output_str

    def save(self, file_path: str) -> None:
        with open(file_path, 'w') as f:
            f.write(self.to_string())

#########################
# Default Config
#########################


LEAP_DEFAULT = {
    "force_field": [
        "protein.ff14SB",
        "gaff2",
        "water.tip3p"
    ],
    "ligand_file_path": None,
    "frcmod_file_path": None,
    "protein_file_path": None,
    "file_prefix": None,
    "box_size": "12.0",
    "box_type": "TIP3PBOX",
    "add_ions_type": "Na+"
}

MINIMIZE_DEFAULT = {
    "imin": 1,
    "cut": 8.0,
    "ntpr": 10,
    "ntb": 1,
    "ntr": 0,
    "ntc": 1,
    "maxcyc": 10000,
    "ncyc": 5000
}

NVT_DEFAULT = {
    "imin": 0,
    "ntx": 5,
    "irest": 1,
    "temp0": 300,
    "ntt": 3,
    "gamma_ln": 2.0,
    "ntp": 0,
    "ntb": 1,
    "ntc": 2,
    "ntf": 2,
    "cut": 8.0,
    "ntpr": 1000,
    "ntwr": 1000,
    "ntwx": 1000,
    "nstlim": 500000,
    "dt": 0.002
}

NPT_DEFAULT = {
    "imin": 0,
    "ntx": 5,
    "irest": 1,
    "temp0": 300,
    "ntt": 3,
    "gamma_ln": 2.0,
    "ntp": 1,
    "taup": 2.0,
    "ntb": 2,
    "ntc": 2,
    "ntf": 2,
    "cut": 8.0,
    "ntpr": 1000,
    "ntwr": 1000,
    "ntwx": 1000,
    "nstlim": 50000000,
    "dt": 0.002
}

#########################
# LEaP Input File Class
#########################


class LeapInput(LeapConstructor):
    def __init__(
            self,
            protein_file_path: str,
            ligand_file_path: Union[str, None] = None,
            frcmod_file_path: Union[str, None] = None,
            file_prefix: Union[str, None] = None,
            box_size: float = 12.0,
            box_type: str = 'TIP3PBOX',
            add_ions_type: str = 'Na+',
            **kwargs) -> None:
        '''
        LEaP Input Config Constructor
        Default Config:
            "force_field": [
                "protein.ff14SB",
                "gaff2",
                "water.tip3p"
                ],
            "ligand_file_path": None,
            "frcmod_file_path": None,
            "protein_file_path": None,
            "file_prefix": None,
            "box_size": "12.0",
            "box_type": "TIP3PBOX",
            "add_ions_type": "Na+"
        '''

        file_prefix = "Untitled" if file_prefix is None else file_prefix
        state_dict = LEAP_DEFAULT.copy()
        state_dict.update(
            {
                "ligand_file_path": ligand_file_path,
                "frcmod_file_path": frcmod_file_path,
                "protein_file_path": protein_file_path,
                "file_prefix": file_prefix,
                "box_size": str(box_size),
                "box_type": box_type,
                "add_ions_type": add_ions_type
            }
        )
        state_dict.update(kwargs)
        super().__init__(state_dict)

#########################
# MD Simulation Input File
#########################


class MinimizeInput(SimulationConstructor):

    def __init__(self, title: str = None, **kwargs) -> None:
        """
        Minimization Input Config Constructor
        Default Attributions:
            {
                "imin": 1,
                "cut": 10.0,
                "ntpr": 10,
                "ntb": 1,
                "ntr": 0,
                "ntc": 1,
                "maxcyc": 10000,
                "ncyc": 5000,
            }
        """
        title = "Minimization" if title is None else title
        state_dict = MINIMIZE_DEFAULT.copy()
        state_dict.update(kwargs)
        super().__init__(state_dict, "cntrl", title=title)


class RestrainedMinimizeInput(SimulationConstructor):
    def __init__(self, restraintmask: str, restraint_wt=2.0, title: str = None, **kwargs) -> None:
        """
        Restrained Minimization Input Config Constructor
        Default Attributions:
            {
                "imin": 1,
                "cut": 10.0,
                "ntpr": 10,
                "ntb": 1,
                "ntr": 1,
                "ntc": 1,
                "maxcyc": 10000,
                "ncyc": 5000,
                "restraintmask": None,
                "restraint_wt": 2.0
            }
        """
        title = 'Restrained Minimization' if title is None else title
        state_dict = MINIMIZE_DEFAULT.copy()
        state_dict.update({
            "ntr": 1,
            "restraintmask": restraintmask,
            "restraint_wt": restraint_wt
        })
        state_dict.update(kwargs)
        super().__init__(state_dict, "cntrl", title=title)


class NVTInput(SimulationConstructor):
    def __init__(self, end_cfg: bool = True, title: str = None, **kwargs) -> None:
        """
        NVT System Input Config Constructor.
        Default Attributions:
            {
                "imin": 0,
                "ntx": 5,
                "irest": 1,
                "temp0": 300,
                "ntt": 3,
                "gamma_ln": 2.0,
                "ntp": 0,
                "ntb": 1,
                "ntc": 2,
                "ntf": 2,
                "cut": 10.0,
                "ntpr": 1000,
                "ntwr": 1000,
                "ntwx": 1000,
                "nstlim": 500000,
                "dt": 0.002
            }

        Parameters
        ----------
        end_cfg: bool, optional
            Whether to add string "END" to the end of input file. Defaults to True.
        title: str, optional
            Title of the simulation. Defaults to None.
        """
        title = "NVT Simulation" if title is None else title
        state_dict = NVT_DEFAULT.copy()
        state_dict.update(kwargs)
        super().__init__(state_dict, "cntrl", end_cfg=end_cfg, title=title)


class HeatInput:
    HEAT_TEMPLATE = {
        "TYPE": "\'TEMP0\'",
        "ISTEP1": 0,
        "ISTEP2": 9000,
        "VALUE1": 0.0,
        "VALUE2": 300.0
    }

    def __init__(
            self, tgt_temperature: float = 300.0, heat_step: int = 9000,
            total_step: int = 10000, step_length: float = 0.002, 
            restraint_wt: float = None, restraintmask: str = None) -> None:
        """
        System Heating Config Constructor.

        Parameters
        ----------
        tgt_temperature: float, optional
            Target temperature of heating. Defaults to 300.0.
        heat_step: int, optional 
            Heating system prior to this step, then, keep temperature until the end. 
            Should not larger than total steps. Defaults to 9000.
        total_step: int, optional
            Total steps of heating progress. Defaults to 10000.
        step_length: float, optional
            Length of each step. Defaults to 0.002.
        restraint_wt: float, optional
            Weight of the restraint. Defaults to None.
        restraintmask: str, optional
            Mask of the restraint. Defaults to None.
        """
        self.stage_1_dict = self.HEAT_TEMPLATE.copy()
        self.stage_1_dict.update({
            "ISTEP2": heat_step,
            "VALUE2": tgt_temperature
        })

        self.stage_2_dict = self.HEAT_TEMPLATE.copy()
        self.stage_2_dict.update({
            "ISTEP1": heat_step + 1,
            "ISTEP2": total_step,
            "VALUE1": tgt_temperature,
            "VALUE2": tgt_temperature
        })
        self.tgt_temperature = tgt_temperature
        self.heat_step = heat_step
        self.total_step = total_step
        self.step_length = step_length
        self.restraint_wt = restraint_wt
        self.restraintmask = restraintmask

        self.manager = MultiConstructorManager()

    @property
    def _state_dict(self) -> dict:
        return self.manager._state_dict
    
    def get_state_dict(self) -> dict:
        return self._state_dict
    
    def add_nvt(self, **kwargs) -> None:
        nvt_part = NVTInput(
            end_cfg=False,
            imin=0, irest=0, ntx=1,
            nstlim=self.total_step, tempi=0.0,
            temp0=self.tgt_temperature,
            dt=self.step_length,
            ntpr=100, ntwr=100, ntwx=100,
            title="Heat",
            **kwargs
        )
        self.manager.add_constructor(nvt_part)

    def add_heat(self, stage_dict: dict = None) -> None:
        if stage_dict is None:
            self.manager.add_constructor(SimulationConstructor(
                self.stage_1_dict, 'wt', end_cfg=False))
            self.manager.add_constructor(SimulationConstructor(
                self.stage_2_dict, 'wt', end_cfg=False))
            self.manager.add_constructor(SimulationConstructor(
                {'TYPE': "\'END\'"}, 'wt', end_cfg=True))
        else:
            self.manager.add_constructor(
                SimulationConstructor(stage_dict, 'wt', end_cfg=False))

    def _default_workflow(self) -> None:
        if self.manager.is_empty():
            if self.restraint_wt is not None and self.restraintmask is not None:
                self.add_nvt(restraintmask=self.restraintmask,
                             restraint_wt=self.restraint_wt)
            else:
                self.add_nvt()
            self.add_heat()

    def to_string(self) -> str:
        self._default_workflow()
        return self.manager.to_string()

    def save(self, file_path: str) -> None:
        self._default_workflow()
        self.manager.save(file_path)


class NPTInput(SimulationConstructor):
    def __init__(self, end_cfg: bool = True, title: str = None, **kwargs) -> None:
        """
        NPT System Input Config Constructor.
        Default Attributions:
            {
                "imin": 0,
                "ntx": 5,
                "irest": 1,
                "temp0": 300,
                "ntt": 3,
                "gamma_ln": 2.0,
                "ntp": 1,
                "taup": 2.0,
                "ntb": 2,
                "ntc": 2,
                "ntf": 2,
                "cut": 10.0,
                "ntpr": 1000,
                "ntwr": 1000,
                "ntwx": 1000,
                "nstlim": 50000000,
                "dt": 0.002
            }
        """
        title = "NPT Simulation" if title is None else title
        state_dict = NPT_DEFAULT.copy()
        state_dict.update(kwargs)
        super().__init__(state_dict, "cntrl", end_cfg=end_cfg, title=title)

#########################
# MD Trajectory Analysis Input File
#########################


MMGBSA_GENERAL_DEFAULT = {
    "startframe": None,
    "endframe": None,
    "step_size": None,
    "verbose": 2,
    "keep_files": 1,
    "netcdf": 1,
}

MMGBSA_GB_DEFAULT = {
    "igb": 5,
    "saltcon": 0.15
}

MMGBSA_DECOMP_DEFAULT = {
    "idecomp": 1,
    "dec_verbose": 1,
}

MMGBSA_PB_DEFAULT = {
    "istrng": 0.15,
    "fillration": 4.0,
}

MMGBSA_NMODE_DEFAULT = {
    "maxcyc": 10000,
    "drms": 0.001,
}


class MMGBSAInput:
    def __init__(
        self,
        start_frame: int, end_frame: int,
        step_size: int = 1
    ) -> None:

        self.general_dict = MMGBSA_GENERAL_DEFAULT.copy()
        self.general_dict.update(
            startframe=start_frame, endframe=end_frame, step_size=step_size)
        self.manager = MultiConstructorManager()

    def add_general(self, title: str = None, **kwargs) -> None:
        self.general_dict.update(kwargs)
        general_part = AnalysisConstructor(
            self.general_dict, 'general', title=title)
        self.manager.add_constructor(general_part)

    def add_gb(self, **kwargs) -> None:
        gb_dict = MMGBSA_GB_DEFAULT.copy()
        gb_dict.update(kwargs)
        gb_part = AnalysisConstructor(gb_dict, 'gb')
        self.manager.add_constructor(gb_part)

    def add_pb(self, **kwargs) -> None:
        pb_dict = MMGBSA_PB_DEFAULT.copy()
        pb_dict.update(kwargs)
        pb_part = AnalysisConstructor(pb_dict, 'pb')
        self.manager.add_constructor(pb_part)

    def add_decomp(self, **kwargs) -> None:
        decomp_dict = MMGBSA_DECOMP_DEFAULT.copy()
        decomp_dict.update(kwargs)
        decomp_part = AnalysisConstructor(decomp_dict, 'decomp')
        self.manager.add_constructor(decomp_part)

    def add_nmode(self, **kwargs) -> None:
        nmode_dict = MMGBSA_NMODE_DEFAULT.copy()
        nmode_dict.update(kwargs)
        nmode_part = AnalysisConstructor(nmode_dict, 'nmode')
        self.manager.add_constructor(nmode_part)

    def to_string(self) -> str:
        return self.manager.to_string()

    def save(self, file_path: str) -> None:
        self.manager.save(file_path)
