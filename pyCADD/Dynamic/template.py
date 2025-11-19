"""Template classes for molecular dynamics input file generation.

This module provides template classes for generating various types of molecular dynamics
input files including Amber MD input files, tLEaP scripts, and MM-GBSA analysis inputs.
All templates use string formatting to replace placeholders with actual simulation parameters.
"""

from abc import ABC, abstractmethod
from typing import Literal

from pyCADD.utils.common import File
from pyCADD.utils.tool import write_file

#########################
# Basic Class Definition
#########################


class BaseConstructor(ABC):
    """Abstract base class for all input file template constructors.

    Provides common functionality for managing template parameters and generating
    formatted input files for molecular dynamics simulations and related tools.

    Attributes:
        _state_dict (dict): Dictionary containing template parameters.
        _state_keys (list): List of parameter keys for state management.
    """

    def __init__(self, state_dict: dict) -> None:
        """Initializes the base constructor with template parameters.

        Args:
            state_dict (dict): Dictionary containing template parameters and values.

        Note:
            All keys in state_dict are automatically set as instance attributes
            for convenient access in template formatting.
        """
        self._state_dict = state_dict.copy()
        self._state_keys = list(self._state_dict.keys())
        for key in self._state_keys:
            self.__dict__[key] = self._state_dict[key]

    def get_state_dict(self) -> dict:
        """Returns the current state of the constructor as a dictionary.

        Returns:
            dict: Dictionary containing current values of all template parameters.

        Note:
            Only returns parameters that were originally in the state_dict,
            not any additional attributes that may have been added later.
        """
        return {key: self.__dict__[key] for key in self._state_keys}

    @abstractmethod
    def to_string(self) -> str:
        """Generates formatted template string with current parameters.

        Returns:
            str: Complete formatted template string ready for file output.

        Note:
            This method must be implemented by all subclasses to define
            how the template should be formatted with the current parameters.
        """
        ...

    def to_str(self) -> str:
        """Alias for to_string() method.

        Returns:
            str: Same as to_string() - formatted template string.
        """
        return self.to_string()

    def to_dict(self) -> dict:
        """Returns the current state as a dictionary.

        Returns:
            dict: Dictionary containing current template parameter values.

        Note:
            This is an alias for get_state_dict() method for convenience.
        """
        return self.get_state_dict()

    def save(self, file_path: str) -> str:
        """Saves the formatted template to a file.

        Args:
            file_path (str): Path where the template file should be saved.

        Returns:
            str: Path to the saved file.

        Note:
            The template is automatically formatted using to_string() before
            saving. Any existing file at the target path will be overwritten.
        """
        return write_file(file_path, self.to_string())

    def add_attr(self, **kwargs) -> None:
        """Adds or updates template attributes.

        Args:
            **kwargs: Attribute name-value pairs to add or update.

        Note:
            Only attributes that were originally in the state dictionary
            can be updated. New attributes not in _state_keys are ignored.
        """
        for key in kwargs.keys():
            if key in self._state_keys:
                raise KeyError(f"Attribute {key} already exists in the current state.")
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
                raise KeyError(f"Attribute {key} can not be found in the current state.")
        self.__dict__.update(kwargs)

    def del_attr(self, attr_name: str) -> None:
        """
        Delete the attribute from the current state.
        """
        if attr_name not in self._state_keys:
            raise KeyError(f"Attribute {attr_name} can not be found in the current state.")
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
{solvatebox} com {box_type} {box_size}
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
    """Template constructor for Amber tLEaP input files.

    Generates tLEaP scripts for building molecular systems including proteins,
    ligands, solvent boxes, and counter ions. Supports both protein-only (apo)
    and protein-ligand systems.

    Required state_dict keys:
        force_field (list): List of force field names (e.g., ["protein.ff14SB", "gaff2"])
        protein_file_path (str): Path to prepared protein PDB file
        file_prefix (str): Prefix for output topology and coordinate files
        box_type (str): Type of solvent box (e.g., "TIP3PBOX")
        box_size (float): Distance from solute to box edge in Angstroms
        add_ions_type (str): Type of neutralizing ions (e.g., "Na+", "Cl-")

    Optional state_dict keys (for protein-ligand systems):
        ligand_file_path (str): Path to ligand MOL2 file
        frcmod_file_path (str): Path to additional force field parameters
    """

    def __init__(self, state_dict: dict) -> None:
        """Initializes the tLEaP template constructor.

        Args:
            state_dict (dict): Dictionary containing template parameters.
                See class docstring for required and optional keys.

        Note:
            Automatically sets pro_lig attribute for complex naming in templates.
            The template selection (apo vs protein-ligand) is determined by the
            presence of ligand_file_path in the state dictionary.
        """
        super().__init__(state_dict)
        self.add_attr(pro_lig=r"{pro lig}")

    @property
    def base_template(self) -> str:
        """Returns the appropriate base template for the system type.

        Returns:
            str: Base tLEaP template string for either apo (protein-only) or
                protein-ligand systems based on the presence of ligand_file_path.
        """
        if self.get_state_dict().get("ligand_file_path") is None:
            return _LEAP_BASE_TEMPLATE_APO
        else:
            return _LEAP_BASE_TEMPLATE

    def _get_force_field_template(self) -> str:
        """Generates force field loading commands for tLEaP.

        Returns:
            str: Multi-line string containing "source leaprc.{ff}" commands
                for all force fields specified in the force_field list.

        Note:
            Force fields are loaded in the order specified in the force_field list.
            Common force fields include protein.ff14SB, gaff2, water.tip3p, etc.
        """
        _template = "\n".join([f"source leaprc.{ff}" for ff in self.get_state_dict()["force_field"]])
        return _template

    def to_string(self) -> str:
        """Generates the complete tLEaP input file content.

        Returns:
            str: Complete tLEaP script with force field loading commands
                followed by the formatted base template.

        Note:
            The output combines force field loading commands with the
            system-specific template (apo or protein-ligand) formatted
            with the current state parameters.
        """
        output_str = self._get_force_field_template()
        output_str += self.base_template.format(**self.get_state_dict())
        return output_str


class SimulationConstructor(BaseConstructor):

    def __init__(
        self,
        state_dict: dict,
        type: Literal["cntrl", "wt"],
        end_cfg: bool = True,
        title: str = None,
    ) -> None:
        """Initialize a SimulationConstructor for MD simulation input files.

        Creates a constructor for generating MD simulation input files used for
        energy minimization, heating, equilibration, and production simulations.
        Supports both control (&cntrl) and weight (&wt) namelist formats.

        Args:
            state_dict: Dictionary containing simulation parameters and their values.
            type: Type of simulation namelist ('cntrl' for standard MD parameters,
                'wt' for weight parameters).
            end_cfg: Whether to add "END" string at the end of input file.
                Defaults to True.
            title: Optional title string for the simulation stage.
                Defaults to None.

        Note:
            The constructor supports two main Amber MD input formats:
            - &cntrl: Standard MD control parameters (imin, nstlim, dt, etc.)
            - &wt: Weight parameters for restraints and temperature ramping
        """
        super().__init__(state_dict)
        self._type = type
        self._end_cfg = end_cfg
        self.title = title

    def to_string(self) -> str:
        """Generate formatted MD simulation input file content.

        Creates a properly formatted Amber MD input file with namelist syntax,
        including optional title, parameter assignments, and termination strings.

        Returns:
            str: Formatted MD input file content with:
                - Optional title line
                - Namelist declaration (&cntrl or &wt)
                - Parameter assignments (key=value format)
                - Namelist termination (/)
                - Optional END statement

        Note:
            Output format follows Amber MD input conventions:
            ```
            [title]
            &namelist_type
            param1=value1,
            param2=value2,
            ...
            /
            [END]
            ```
        """
        if self.title is not None:
            output_str = f"{self.title}\n&{self._type}\n"
        else:
            output_str = f"&{self._type}\n"
        output_str += ",\n".join(
            [f"{str(key)}={str(value)}" for key, value in self.get_state_dict().items()]
        )
        output_str += "\n/\n"
        if self._end_cfg:
            output_str += "\nEND"
        return output_str


class AnalysisConstructor(BaseConstructor):
    """Constructor for MD trajectory analysis input files.

    Generates input files for various types of MD trajectory analysis including
    general analysis, GB/PB calculations, decomposition analysis, and normal
    mode analysis using Amber analysis tools.

    Attributes:
        _type: Type of analysis to perform.
        title: Optional title for the analysis.
    """

    def __init__(
        self,
        state_dict: dict,
        type: Literal["general", "gb", "pb", "decomp", "nmode"],
        title: str = None,
    ) -> None:
        """Initialize an AnalysisConstructor for MD trajectory analysis.

        Creates a constructor for generating analysis input files used with
        Amber trajectory analysis tools like cpptraj, MMPBSA, and nmode.

        Args:
            state_dict: Dictionary containing analysis parameters and their values.
            type: Type of analysis to perform:
                - 'general': General trajectory analysis (RMSD, RMSF, etc.)
                - 'gb': Generalized Born solvation analysis
                - 'pb': Poisson-Boltzmann solvation analysis  
                - 'decomp': Energy decomposition analysis
                - 'nmode': Normal mode analysis
            title: Optional title string for the analysis job.
                Defaults to None.

        Note:
            Different analysis types require specific parameter sets and
            are processed by different Amber analysis programs.
        """
        super().__init__(state_dict)
        self._type = type
        self.title = title

    def to_string(self) -> str:
        """Generate formatted analysis input file content.

        Creates a properly formatted Amber analysis input file with comment
        headers, namelist syntax, and parameter assignments for trajectory
        analysis tools.

        Returns:
            str: Formatted analysis input file content with:
                - Optional comment title line (# title)
                - Namelist declaration (&analysis_type)
                - Parameter assignments with leading space
                - Namelist termination (/)

        Note:
            Output format follows Amber analysis tool conventions:
            ```
            # [title]
            &analysis_type
             param1=value1,
             param2=value2,
            ...
            /
            ```
        """
        output_str = f"# {self.title}\n" if self.title is not None else ""
        output_str += f"&{self._type}\n"
        output_str += ",\n".join(
            [f" {str(key)}={str(value)}" for key, value in self.get_state_dict().items()]
        )
        output_str += "\n/\n"
        return output_str


class MultiConstructorManager:
    """Manager for handling multiple template constructors.

    Provides a unified interface for managing multiple BaseConstructor instances,
    allowing batch operations, state aggregation, and combined output generation
    for complex MD workflows requiring multiple input files.

    Attributes:
        constructor_list: List of BaseConstructor instances managed by this manager.
    """

    def __init__(self, constructor_list: list = None) -> None:
        """Initialize a MultiConstructorManager.

        Args:
            constructor_list: Optional list of BaseConstructor instances to manage.
                Defaults to empty list if None.
        """
        self.constructor_list = constructor_list if constructor_list is not None else []

    @property
    def _state_dict(self) -> dict:
        """Get aggregated state dictionary from all managed constructors.

        Combines state dictionaries from all constructors in the manager,
        with later constructors potentially overriding parameters from
        earlier ones.

        Returns:
            dict: Aggregated state dictionary containing all parameters
                from managed constructors. Empty dict if no constructors present.

        Note:
            Parameter conflicts are resolved by last-constructor-wins policy.
        """
        _state_dict = dict()
        if len(self.constructor_list) == 0:
            return _state_dict
        else:
            for constructor in self.constructor_list:
                _state_dict.update(constructor.get_state_dict())
            return _state_dict

    def is_empty(self) -> bool:
        """Check if the manager contains any constructors.

        Returns:
            bool: True if constructor list is empty, False otherwise.
        """
        return len(self.constructor_list) == 0

    def add_constructor(self, constructor: BaseConstructor) -> None:
        """Add a constructor to the management list.

        Args:
            constructor: BaseConstructor instance to add to the manager.
        """
        self.constructor_list.append(constructor)

    def del_constructor(self, constructor: BaseConstructor) -> None:
        """Remove a constructor from the management list.

        Args:
            constructor: BaseConstructor instance to remove from the manager.

        Raises:
            ValueError: If constructor is not found in the list.
        """
        self.constructor_list.remove(constructor)

    def to_string(self) -> str:
        """Generate combined string output from all managed constructors.

        Concatenates the string output from all constructors in the management
        list, maintaining the order they were added.

        Returns:
            str: Combined string output from all managed constructors.
                Empty string if no constructors present.
        """
        output_str = "".join([tmp.to_string() for tmp in self.constructor_list])
        return output_str

    def save(self, file_path: str) -> str:
        """Save combined constructor output to file.

        Writes the concatenated string output from all managed constructors
        to the specified file path.

        Args:
            file_path: Path where the combined output should be saved.

        Returns:
            str: Path to the saved file.

        Note:
            Uses the write_file utility function for file operations.
        """
        return write_file(file_path, self.to_string())


#########################
# Default Config
#########################


LEAP_DEFAULT = {
    "force_field": ["protein.ff14SB", "gaff2", "water.tip3p"],
    "ligand_file_path": None,
    "frcmod_file_path": None,
    "protein_file_path": None,
    "file_prefix": None,
    "box_size": "12.0",
    "box_type": "TIP3PBOX",
    "add_ions_type": "Na+",
    "solvatebox": "solvatebox"
}

MINIMIZE_DEFAULT = {
    "imin": 1,
    "cut": 8.0,
    "ntpr": 10,
    "ntb": 1,
    "ntr": 0,
    "ntc": 1,
    "maxcyc": 10000,
    "ncyc": 5000,
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
    "dt": 0.002,
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
    "dt": 0.002,
}

#########################
# LEaP Input File Class
#########################


class LeapInput(LeapConstructor):
    """tLEaP input file constructor for system preparation.

    Specialized constructor for generating tLEaP scripts that prepare MD systems
    including protein structures, ligands, solvation, and ion neutralization.
    Extends LeapConstructor with convenient parameter management for common
    MD system setup tasks.
    """

    def __init__(
        self,
        protein_file_path: str | File,
        ligand_file_path: str | File = None,
        frcmod_file_path: str | File = None,
        file_prefix: str = None,
        box_size: float = 10.0,
        box_type: str = "TIP3PBOX",
        solvatebox: str = "solvatebox",
        add_ions_type: str = "Na+",
        **kwargs,
    ) -> None:
        """Initialize a LeapInput constructor for tLEaP system preparation.

        Creates a tLEaP input constructor with default force fields and
        customizable system parameters for protein-ligand MD simulations.

        Args:
            protein_file_path: Path to protein structure file (PDB format).
            ligand_file_path: Optional path to ligand structure file.
                Defaults to None for apo simulations.
            frcmod_file_path: Optional path to custom force field modification file.
                Defaults to None.
            file_prefix: Prefix for generated output files.
                Defaults to "Untitled" if None.
            box_size: Size of solvation box in Angstroms.
                Defaults to 10.0.
            box_type: Type of water box (TIP3PBOX, TIP4PBOX, etc.).
                Defaults to "TIP3PBOX".
            solvatebox: Solvation command type.
                Defaults to "solvatebox".
            add_ions_type: Type of ions for neutralization.
                Defaults to "Na+".
            **kwargs: Additional parameters to override default values.

        Note:
            Default force field combination includes:
            - protein.ff14SB: Amber protein force field
            - gaff2: General Amber force field for small molecules  
            - water.tip3p: TIP3P water model
        """

        file_prefix = "Untitled" if file_prefix is None else file_prefix
        state_dict = LEAP_DEFAULT.copy()
        state_dict.update(
            {
                "ligand_file_path": File(ligand_file_path).file_path if ligand_file_path else None,
                "frcmod_file_path": File(frcmod_file_path).file_path if frcmod_file_path else None,
                "protein_file_path": File(protein_file_path).file_path,
                "file_prefix": file_prefix,
                "box_size": str(box_size),
                "box_type": box_type,
                "add_ions_type": add_ions_type,
                "solvatebox": solvatebox,
            }
        )
        state_dict.update(kwargs)
        super().__init__(state_dict)


#########################
# MD Simulation Input File
#########################


class MinimizeInput(SimulationConstructor):
    """Constructor for energy minimization input files.

    Generates Amber input files for energy minimization simulations using
    steepest descent and conjugate gradient algorithms to relax molecular
    structures before MD simulations.
    """

    def __init__(self, title: str = None, **kwargs) -> None:
        """Initialize a MinimizeInput constructor.

        Creates an energy minimization input with default parameters optimized
        for typical protein-ligand systems.

        Args:
            title: Optional title for the minimization job.
                Defaults to "Minimization".
            **kwargs: Additional parameters to override defaults.

        Note:
            Default minimization parameters:
            - imin=1: Energy minimization flag
            - cut=8.0: Nonbonded cutoff distance (Å)
            - ntpr=10: Print frequency
            - ntb=1: Constant volume periodic boundaries
            - ntr=0: No positional restraints
            - ntc=1: No SHAKE constraints
            - maxcyc=10000: Maximum minimization cycles
            - ncyc=5000: Steepest descent cycles before conjugate gradient
        """
        title = "Minimization" if title is None else title
        state_dict = MINIMIZE_DEFAULT.copy()
        state_dict.update(kwargs)
        super().__init__(state_dict, "cntrl", title=title)


class RestrainedMinimizeInput(SimulationConstructor):
    """Constructor for restrained energy minimization input files.

    Generates Amber input files for energy minimization with positional
    restraints on specified atoms, typically used to minimize solvent
    and ions while keeping protein structure fixed.
    """

    def __init__(self, restraintmask: str, restraint_wt=2.0, title: str = None, **kwargs) -> None:
        """Initialize a RestrainedMinimizeInput constructor.

        Creates a restrained energy minimization input with specified atom
        restraints and force constants.

        Args:
            restraintmask: Amber mask syntax specifying atoms to restrain
                (e.g., ":1-100@CA,C,N,O" for protein backbone).
            restraint_wt: Force constant for positional restraints in kcal/mol/Å².
                Defaults to 2.0.
            title: Optional title for the minimization job.
                Defaults to "Restrained Minimization".
            **kwargs: Additional parameters to override defaults.

        Note:
            Builds on standard minimization with:
            - ntr=1: Enable positional restraints
            - restraintmask: Atom selection for restraints
            - restraint_wt: Restraint force constant
        """
        title = "Restrained Minimization" if title is None else title
        state_dict = MINIMIZE_DEFAULT.copy()
        state_dict.update({"ntr": 1, "restraintmask": restraintmask, "restraint_wt": restraint_wt})
        state_dict.update(kwargs)
        super().__init__(state_dict, "cntrl", title=title)


class NVTInput(SimulationConstructor):
    """Constructor for NVT (constant volume, temperature) MD simulation input files.

    Generates Amber input files for NVT ensemble molecular dynamics simulations
    with temperature control via Langevin dynamics, typically used for
    equilibration after minimization or heating phases.
    """

    def __init__(self, end_cfg: bool = True, title: str = None, **kwargs) -> None:
        """Initialize an NVTInput constructor.

        Creates an NVT simulation input with default parameters for temperature
        equilibration at 300 K using Langevin thermostat.

        Args:
            end_cfg: Whether to add "END" string at end of input file.
                Defaults to True.
            title: Optional title for the simulation.
                Defaults to "NVT Simulation".
            **kwargs: Additional parameters to override defaults.

        Note:
            Default NVT parameters:
            - imin=0: MD simulation flag
            - ntx=5: Read coordinates and velocities
            - irest=1: Restart simulation
            - temp0=300: Target temperature (K)
            - ntt=3: Langevin thermostat
            - gamma_ln=2.0: Collision frequency (ps⁻¹)
            - ntp=0: No pressure control
            - ntb=1: Constant volume periodic boundaries
            - ntc=2: SHAKE constraints on H-bonds
            - ntf=2: Omit forces for constrained bonds
            - nstlim=500000: Simulation steps (1 ns at dt=0.002 ps)
        """
        title = "NVT Simulation" if title is None else title
        state_dict = NVT_DEFAULT.copy()
        state_dict.update(kwargs)
        super().__init__(state_dict, "cntrl", end_cfg=end_cfg, title=title)


class HeatInput:
    """Constructor for system heating MD simulation with temperature ramping.

    Manages gradual temperature increase from 0 K to target temperature using
    Amber weight (&wt) namelist for controlled heating phases. Combines NVT
    simulation parameters with temperature ramping schedules.

    Attributes:
        HEAT_TEMPLATE: Default weight parameters for temperature ramping.
        stage_1_dict: Weight parameters for heating phase.
        stage_2_dict: Weight parameters for temperature maintenance phase.
        manager: MultiConstructorManager for combining simulation components.
    """

    HEAT_TEMPLATE = {"TYPE": "'TEMP0'", "ISTEP1": 0, "ISTEP2": 9000, "VALUE1": 0.0, "VALUE2": 300.0}

    def __init__(
        self,
        target_temp: float = 300.0,
        heat_step: int = 9000,
        total_step: int = 10000,
        step_length: float = 0.002,
        restraint_wt: float = None,
        restraintmask: str = None,
    ) -> None:
        """Initialize a HeatInput constructor for controlled heating.

        Creates a heating protocol with linear temperature ramping followed
        by temperature maintenance using Amber weight parameters.

        Args:
            target_temp: Target temperature for heating in Kelvin.
                Defaults to 300.0.
            heat_step: Number of steps for temperature ramping phase.
                Must be less than total_step. Defaults to 9000.
            total_step: Total simulation steps for entire heating protocol.
                Defaults to 10000.
            step_length: Time step length in picoseconds.
                Defaults to 0.002.
            restraint_wt: Optional restraint force constant in kcal/mol/Å².
                Defaults to None.
            restraintmask: Optional Amber mask for restrained atoms.
                Defaults to None.

        Note:
            Creates two-stage heating protocol:
            1. Ramping: 0 K → target_temp over heat_step steps
            2. Maintenance: target_temp held for remaining steps
        """
        self.stage_1_dict = self.HEAT_TEMPLATE.copy()
        self.stage_1_dict.update({"ISTEP2": heat_step, "VALUE2": target_temp})

        self.stage_2_dict = self.HEAT_TEMPLATE.copy()
        self.stage_2_dict.update(
            {
                "ISTEP1": heat_step + 1,
                "ISTEP2": total_step,
                "VALUE1": target_temp,
                "VALUE2": target_temp,
            }
        )
        self.target_temp = target_temp
        self.heat_step = heat_step
        self.total_step = total_step
        self.step_length = step_length
        self.restraint_wt = restraint_wt
        self.restraintmask = restraintmask

        self.manager = MultiConstructorManager()

    @property
    def _state_dict(self) -> dict:
        """Get aggregated state dictionary from managed constructors.

        Returns:
            dict: Combined state parameters from all managed constructors.
        """
        return self.manager._state_dict

    def get_state_dict(self) -> dict:
        """Get current state dictionary.

        Returns:
            dict: Current state parameters for the heating protocol.
        """
        return self._state_dict

    def add_nvt(self, **kwargs) -> None:
        """Add NVT simulation component to heating protocol.

        Creates and adds an NVT simulation constructor configured for heating
        with appropriate initial conditions and output frequencies.

        Args:
            **kwargs: Additional parameters to override NVT defaults.

        Note:
            Configures NVT for heating with:
            - irest=0, ntx=1: Cold start from coordinates only
            - tempi=0.0: Initial temperature 0 K
            - Custom output frequencies for heating monitoring
        """
        nvt_part = NVTInput(
            end_cfg=False,
            imin=0,
            irest=0,
            ntx=1,
            nstlim=self.total_step,
            tempi=0.0,
            temp0=self.target_temp,
            dt=self.step_length,
            ntpr=100,
            ntwr=100,
            ntwx=100,
            title="Heat",
            **kwargs,
        )
        self.manager.add_constructor(nvt_part)

    def add_heat(self, stage_dict: dict = None) -> None:
        """Add temperature ramping weight parameters to heating protocol.

        Adds weight (&wt) namelist sections for controlling temperature ramping
        during the heating simulation.

        Args:
            stage_dict: Optional custom weight parameters dictionary.
                If None, uses default two-stage heating protocol.

        Note:
            Default protocol adds:
            1. Stage 1: Linear ramp from 0 K to target temperature
            2. Stage 2: Hold at target temperature
            3. END directive to terminate weight section
        """
        if stage_dict is None:
            self.manager.add_constructor(
                SimulationConstructor(self.stage_1_dict, "wt", end_cfg=False)
            )
            self.manager.add_constructor(
                SimulationConstructor(self.stage_2_dict, "wt", end_cfg=False)
            )
            self.manager.add_constructor(
                SimulationConstructor({"TYPE": "'END'"}, "wt", end_cfg=True)
            )
        else:
            self.manager.add_constructor(SimulationConstructor(stage_dict, "wt", end_cfg=False))

    def _default_workflow(self) -> None:
        """Execute default heating workflow if manager is empty.

        Automatically sets up complete heating protocol with NVT simulation
        and temperature ramping if no constructors have been added manually.

        Note:
            Applies restraints if both restraint_wt and restraintmask are provided.
        """
        if self.manager.is_empty():
            if self.restraint_wt is not None and self.restraintmask is not None:
                self.add_nvt(restraintmask=self.restraintmask, restraint_wt=self.restraint_wt)
            else:
                self.add_nvt()
            self.add_heat()

    def to_string(self) -> str:
        """Generate complete heating input file content.

        Executes default workflow if needed, then returns formatted input
        file content for Amber heating simulation.

        Returns:
            str: Complete Amber input file with NVT and weight sections.
        """
        self._default_workflow()
        return self.manager.to_string()

    def save(self, file_path: str) -> None:
        """Save heating input file to disk.

        Executes default workflow if needed, then saves complete input
        file content to the specified path.

        Args:
            file_path: Path where the heating input file should be saved.
        """
        self._default_workflow()
        self.manager.save(file_path)


class NPTInput(SimulationConstructor):
    """Constructor for NPT (constant pressure, temperature) MD simulation input files.

    Generates Amber input files for NPT ensemble molecular dynamics simulations
    with both temperature and pressure control, typically used for production
    simulations and system equilibration at physiological conditions.
    """

    def __init__(self, end_cfg: bool = True, title: str = None, **kwargs) -> None:
        """Initialize an NPTInput constructor.

        Creates an NPT simulation input with default parameters for production
        MD simulations at 300 K and 1 atm using Langevin thermostat and
        Berendsen barostat.

        Args:
            end_cfg: Whether to add "END" string at end of input file.
                Defaults to True.
            title: Optional title for the simulation.
                Defaults to "NPT Simulation".
            **kwargs: Additional parameters to override defaults.

        Note:
            Default NPT parameters:
            - imin=0: MD simulation flag
            - ntx=5: Read coordinates and velocities
            - irest=1: Restart simulation
            - temp0=300: Target temperature (K)
            - ntt=3: Langevin thermostat
            - ntp=1: Isotropic pressure scaling
            - taup=2.0: Pressure relaxation time (ps)
            - ntb=2: Constant pressure periodic boundaries
            - nstlim=50000000: Long production simulation (100 ns)
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
    "interval": None,
    "verbose": 2,
    "keep_files": 1,
    "netcdf": 1,
}

MMGBSA_GB_DEFAULT = {"igb": 5, "saltcon": 0.15}

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
    """Constructor for MM/GBSA (MM/PBSA) analysis input files.

    Generates input files for Molecular Mechanics/Generalized Born Surface Area
    or Poisson-Boltzmann Surface Area free energy calculations using AMBER's
    MMPBSA.py tool. Supports various analysis methods including GB, PB,
    decomposition, and normal mode analysis.

    Attributes:
        general_dict: General analysis parameters.
        manager: MultiConstructorManager for combining analysis components.
    """

    def __init__(self, start_frame: int, end_frame: int, step_size: int = 1) -> None:
        """Initialize an MMGBSAInput constructor.

        Sets up MM/GBSA analysis with frame range specification and default
        general parameters for trajectory processing.

        Args:
            start_frame: First trajectory frame to analyze (1-indexed).
            end_frame: Last trajectory frame to analyze (inclusive).
            step_size: Frame interval for analysis sampling.
                Defaults to 1 (analyze every frame).

        Note:
            Frame numbering follows AMBER convention (1-indexed).
            Analysis will process frames: start_frame, start_frame+step_size,
            ..., up to end_frame.
        """

        self.general_dict = MMGBSA_GENERAL_DEFAULT.copy()
        self.general_dict.update(startframe=start_frame, endframe=end_frame, interval=step_size)
        self.manager = MultiConstructorManager()

    def add_general(self, title: str = None, **kwargs) -> None:
        """Add general analysis parameters section.

        Adds the &general namelist with trajectory processing parameters,
        file handling options, and output control settings.

        Args:
            title: Optional title comment for the analysis.
            **kwargs: Additional general parameters to override defaults.

        Note:
            Common parameters include verbose, keep_files, netcdf for
            controlling output detail and file management.
        """
        self.general_dict.update(kwargs)
        general_part = AnalysisConstructor(self.general_dict, "general", title=title)
        self.manager.add_constructor(general_part)

    def add_gb(self, **kwargs) -> None:
        """Add Generalized Born analysis parameters section.

        Adds the &gb namelist for MM/GBSA calculations using implicit
        solvent models based on Generalized Born theory.

        Args:
            **kwargs: GB parameters to override defaults.

        Note:
            Default GB model (igb=5) with physiological salt concentration
            (saltcon=0.15 M). Other igb values: 1, 2, 5, 7, 8 for different
            GB models (Hawkins, Onufriev, etc.).
        """
        gb_dict = MMGBSA_GB_DEFAULT.copy()
        gb_dict.update(kwargs)
        gb_part = AnalysisConstructor(gb_dict, "gb")
        self.manager.add_constructor(gb_part)

    def add_pb(self, **kwargs) -> None:
        """Add Poisson-Boltzmann analysis parameters section.

        Adds the &pb namelist for MM/PBSA calculations using explicit
        solution of the Poisson-Boltzmann equation for electrostatics.

        Args:
            **kwargs: PB parameters to override defaults.

        Note:
            Default parameters include ionic strength (istrng=0.15) and
            grid fill ratio (fillration=4.0) for finite difference
            Poisson-Boltzmann solver.
        """
        pb_dict = MMGBSA_PB_DEFAULT.copy()
        pb_dict.update(kwargs)
        pb_part = AnalysisConstructor(pb_dict, "pb")
        self.manager.add_constructor(pb_part)

    def add_decomp(self, **kwargs) -> None:
        """Add energy decomposition analysis parameters section.

        Adds the &decomp namelist for per-residue energy decomposition
        analysis to identify key binding site residues and interactions.

        Args:
            **kwargs: Decomposition parameters to override defaults.

        Note:
            Default enables residue-based decomposition (idecomp=1) with
            verbose output (dec_verbose=1) for detailed interaction analysis.
        """
        decomp_dict = MMGBSA_DECOMP_DEFAULT.copy()
        decomp_dict.update(kwargs)
        decomp_part = AnalysisConstructor(decomp_dict, "decomp")
        self.manager.add_constructor(decomp_part)

    def add_nmode(self, **kwargs) -> None:
        """Add normal mode analysis parameters section.

        Adds the &nmode namelist for calculating conformational entropy
        contributions using normal mode analysis of energy-minimized
        structures.

        Args:
            **kwargs: Normal mode parameters to override defaults.

        Note:
            Default parameters include convergence criteria (maxcyc=10000,
            drms=0.001) for minimization prior to normal mode calculation.
            Computationally expensive - use sparingly.
        """
        nmode_dict = MMGBSA_NMODE_DEFAULT.copy()
        nmode_dict.update(kwargs)
        nmode_part = AnalysisConstructor(nmode_dict, "nmode")
        self.manager.add_constructor(nmode_part)

    def to_string(self) -> str:
        """Generate complete MM/GBSA input file content.

        Returns:
            str: Formatted MMPBSA.py input file with all added analysis sections.
        """
        return self.manager.to_string()

    def save(self, file_path: str) -> None:
        """Save MM/GBSA input file to disk.

        Args:
            file_path: Path where the MMPBSA.py input file should be saved.
        """
        self.manager.save(file_path)
