"""Core molecular dynamics simulation preparation and execution functions.

This module provides the main functionality for preparing molecular systems and
running molecular dynamics simulations using Amber tools. It includes protein
preparation, ligand parameterization, system building, and simulation execution.
"""

import logging
import os
import re
from typing import Literal

from pyCADD.Dynamic.constant import CPU_NUM, DEFAULT_STRIP_MASKS, PMEMD, SANDER
from pyCADD.Dynamic.template import LeapInput, MMGBSAInput
from pyCADD.Dynamic.utils import _trace_progress, calc_am1bcc, run_parmchk2
from pyCADD.utils.common import ChDir, File
from pyCADD.utils.tool import is_mpirun_available, read_file, shell_run, timeit, write_file

logger = logging.getLogger(__name__)


def protein_prepare(
    protein_file: File | str,
    save_dir: str = None,
    keep_water: bool = False,
    overwrite: bool = False,
) -> File:
    """Prepares protein structure for molecular dynamics simulation.

    Processes a protein PDB file using pdb4amber to clean up the structure,
    add missing atoms, and prepare it for Amber force field parameterization.
    The preparation includes removing/adding atoms and fixing structural issues.

    Args:
        protein_file (File | str): Input protein PDB file to be processed.
        save_dir (str, optional): Directory to save processed files.
            Defaults to current working directory.
        keep_water (bool, optional): Whether to retain water molecules in final structure.
            Defaults to False.
        overwrite (bool, optional): Whether to overwrite existing prepared files.
            Defaults to False.

    Returns:
        File: Processed protein file ready for MD simulation setup.

    Note:
        The preparation process involves multiple steps:
        1. Remove water and add missing atoms (dry step)
        2. Remove all hydrogens (noH step)
        3. Final cleanup and preparation
        If keep_water=True, water molecules from original file are preserved.

    Raises:
        FileNotFoundError: If the input protein file doesn't exist.
        RuntimeError: If pdb4amber processing fails.
    """
    protein_file = File(protein_file)
    save_dir = save_dir or os.getcwd()
    os.makedirs(save_dir, exist_ok=True)
    file_path = protein_file.file_path
    file_prefix = protein_file.file_prefix

    dry_file_path = os.path.join(save_dir, file_prefix + "_dry.pdb")
    noH_file_path = os.path.join(save_dir, file_prefix + "_noH.pdb")
    prepared_file_path = os.path.join(save_dir, file_prefix + "_prepared.pdb")
    if not overwrite and os.path.exists(prepared_file_path):
        return File(prepared_file_path)

    shell_run(f"pdb4amber -i {file_path} -o {dry_file_path} -p --dry --add-missing-atoms ")
    shell_run(f"pdb4amber -i {dry_file_path} -o {noH_file_path} -y")
    shell_run(f"pdb4amber -i {noH_file_path} -o {prepared_file_path}")

    if keep_water:
        with ChDir(path=save_dir):
            shell_run(f"pdb4amber -i {file_path} -o tmp.pdb")
            shell_run(f'cat tmp.pdb | grep "HOH" > water.pdb')
            shell_run(
                f'cat {prepared_file_path} | grep -v "END" > tmp.pdb && \
                cat water.pdb | sed "s/HOH/WAT/" >> tmp.pdb && \
                echo END >> tmp.pdb'
            )
            shell_run(f"pdb4amber -i tmp.pdb -o {prepared_file_path}")
            shell_run(f"rm tmp.pdb water.pdb")
    return File(prepared_file_path)


def molecule_prepare(
    ligand_file: File,
    charge_method: Literal["resp", "bcc", "resp2"] = "resp",
    charge: int = None,
    multiplicity: int = None,
    dft: str = "B3LYP",
    basis_set: str = "6-31g*",
    cpu_num: int = 4,
    mem_use: str = "8GB",
    solvent: str = "water",
    delta: float = 0.5,
    save_dir: str = None,
    overwrite: bool = False,
) -> tuple[File, File]:
    """Prepares small molecule/ligand for molecular dynamics simulation.

    Generates Amber-compatible parameter and coordinate files for a small molecule
    by calculating partial charges using quantum mechanical methods (RESP, AM1-BCC,
    or RESP2) and creating force field parameters using GAFF.

    Args:
        ligand_file (File): Input ligand structure file (SDF, MOL2, PDB, etc.).
        charge_method (Literal["resp", "bcc", "resp2"], optional): Method for charge calculation.
            Defaults to "resp".
        charge (int, optional): Net charge of the molecule. If None, will be determined automatically.
        multiplicity (int, optional): Spin multiplicity. If None, will be determined automatically.
        dft (str, optional): DFT functional for quantum calculations. Defaults to "B3LYP".
        basis_set (str, optional): Basis set for quantum calculations. Defaults to "6-31g*".
        cpu_num (int, optional): Number of CPUs for quantum calculations. Defaults to 4.
        mem_use (str, optional): Memory allocation for quantum calculations. Defaults to "8GB".
        solvent (str, optional): Solvent for implicit solvation in quantum calculations.
            Defaults to "water".
        delta (float, optional): Interpolation parameter for RESP2 method. Defaults to 0.5.
        save_dir (str, optional): Directory to save generated files.
            Defaults to current working directory.
        overwrite (bool, optional): Whether to overwrite existing parameter files.
            Defaults to False.

    Returns:
        tuple[File, File]: Tuple containing (parameter_file, coordinate_file) where:
            - parameter_file: Amber parameter file (.prmtop)
            - coordinate_file: Amber coordinate file (.inpcrd)

    Raises:
        ValueError: If unsupported charge method is specified.
        RuntimeError: If antechamber or other Amber tools fail.
        FileNotFoundError: If input ligand file doesn't exist.

    Note:
        The function automatically handles:
        - Charge and multiplicity determination if not provided
        - GAFF parameter assignment using antechamber
        - Missing parameter generation using parmchk2
        - System building using tleap
    """

    save_dir = save_dir or os.getcwd()
    os.makedirs(save_dir, exist_ok=True)
    file_path = ligand_file.file_path
    file_prefix = f"{ligand_file.file_prefix}_{charge_method.lower()}"
    cpu_num = cpu_num or CPU_NUM
    charge = charge or 0
    multiplicity = multiplicity or 1
    solvent = solvent or "water"
    mol2_file_path = os.path.join(save_dir, f"{file_prefix}.mol2")

    if os.path.exists(mol2_file_path) and not overwrite:
        return File(mol2_file_path), run_parmchk2(mol2_file_path, save_dir=save_dir)

    logger.info("Calculating the partial charge of ligand...")
    if charge_method.lower() == "bcc":
        calc_method = calc_am1bcc
        kwargs = {
            "save_dir": save_dir,
            "charge": charge,
            "multiplicity": multiplicity,
        }
    else:
        from pyCADD.Density.base import Gauss

        gauss = Gauss()
        kwargs = {
            "charge": charge,
            "multiplicity": multiplicity,
            "dft": dft,
            "basis_set": basis_set,
            "mem_use": mem_use,
            "cpu_num": cpu_num,
            "save_dir": save_dir,
        }
        if charge_method.lower() == "resp2":
            calc_method = gauss.calc_resp2
            kwargs.update({"delta": delta})
        elif charge_method.lower() == "resp":
            calc_method = gauss.calc_resp
            kwargs.update({"solvent": solvent})
        else:
            raise ValueError(f"Unsupported charge calculation type: {charge_method}")
    charged_mol2_file = File(calc_method(file_path, **kwargs))
    frcmod_file = run_parmchk2(charged_mol2_file, save_dir=save_dir)
    return charged_mol2_file, frcmod_file


def _create_leap_inputfile(
    protein_file_path: str,
    ligand_file_path: File | str = None,
    frcmod_file_path: File | str = None,
    file_prefix: str = None,
    box_size: float = 10.0,
    box_type: str = "TIP3PBOX",
    solvatebox: str = "solvatebox",
    save_dir: str = None,
) -> File:
    """Creates a tLEaP input file for system preparation.

    Generates a tLEaP input script for building an Amber molecular system
    including protein, ligand (optional), solvent box, and counter ions.

    Args:
        protein_file_path (str): Path to the prepared protein PDB file.
        ligand_file_path (File | str, optional): Path to the ligand MOL2 file.
            Defaults to None for protein-only systems.
        frcmod_file_path (File | str, optional): Path to additional force field
            parameters file. Defaults to None.
        file_prefix (str, optional): Prefix for output files. Defaults to "Untitled".
        box_size (float, optional): Distance from solute to box edge in Angstroms.
            Defaults to 10.0.
        box_type (str, optional): Type of water box model. Defaults to "TIP3PBOX".
        solvatebox (str, optional): Solvation command type. Defaults to "solvatebox".
        save_dir (str, optional): Directory to save the input file.
            Defaults to current working directory.

    Returns:
        File: Generated tLEaP input file object.

    Note:
        The generated input file includes commands for:
        - Loading appropriate force fields
        - Loading protein and ligand structures
        - Adding solvent box and neutralizing ions
        - Saving parameter and coordinate files
    """
    file_prefix = file_prefix or "Untitled"
    protein_file = File(protein_file_path)
    ligand_file = File(ligand_file_path) if ligand_file_path else None
    frcmod_file = File(frcmod_file_path) if frcmod_file_path else None
    save_dir = save_dir or os.getcwd()
    os.makedirs(save_dir, exist_ok=True)
    input_file_path = os.path.join(save_dir, f"{file_prefix}_tleap.in")
    leap_input = LeapInput(
        protein_file_path=protein_file,
        ligand_file_path=ligand_file,
        frcmod_file_path=frcmod_file,
        file_prefix=file_prefix,
        box_size=box_size,
        box_type=box_type,
        solvatebox=solvatebox,
    )
    leap_input.save(input_file_path)
    return File(input_file_path)


def _get_strip_prmtop_str(prmtop_file: File, strip_masks: str, save_dir: str = None) -> str:
    with ChDir(path=save_dir):
        stripped_file_path = f"{prmtop_file.file_prefix}_stripped.prmtop"
        strip_cmd = f"parmed {prmtop_file.file_path} << EOF\n"
        strip_cmd += f"strip {strip_masks}\n"
        strip_cmd += f"parmout {stripped_file_path}\nEOF"
        shell_run(strip_cmd)
        return File(stripped_file_path).read()


def leap_prepare(
    protein_file: File,
    ligand_file: File = None,
    frcmod_file: File = None,
    box_size: float = 10.0,
    box_type: str = "TIP3PBOX",
    solvatebox: str = "solvatebox",
    save_dir: str = None,
) -> tuple[File, File, File, File, File, File | None]:
    """Prepares solvated molecular system using tLEaP.

    Builds a complete molecular dynamics system by combining protein, ligand (optional),
    solvent, and counter ions using Amber's tLEaP program. Generates topology and
    coordinate files ready for MD simulation.

    Args:
        protein_file (File): Prepared protein structure file.
        ligand_file (File, optional): Ligand MOL2 file with parameters. Defaults to None.
        frcmod_file (File, optional): Additional force field parameter file. Defaults to None.
        box_size (float, optional): Distance from solute to box boundary in Angstroms.
            Defaults to 10.0.
        box_type (str, optional): Water model for solvation. Defaults to "TIP3PBOX".
        solvatebox (str, optional): Solvation command type. Defaults to "solvatebox".
        save_dir (str, optional): Directory to save output files.
            Defaults to current working directory.

    Returns:
        tuple: Tuple containing (topology_file, coordinate_file, pdb_file, com_topfile, pro_topfile, lig_topfile) where:
            - topology_file: Amber parameter/topology file (.prmtop)
            - coordinate_file: Amber coordinate file (.inpcrd)
            - pdb_file: PDB file of the solvated system
            - com_topfile: Combined system topology file (.prmtop)
            - pro_topfile: Protein-only topology file (.prmtop)
            - lig_topfile: Ligand-only topology file (.prmtop) or None if no ligand

    Raises:
        RuntimeError: If tLEaP execution fails.
        FileNotFoundError: If input files don't exist.

    Note:
        The function automatically:
        - Generates appropriate file prefixes based on input files
        - Creates tLEaP input script and executes it
        - Converts coordinates to PDB format for visualization
        - Handles both protein-only and protein-ligand systems
    """
    save_dir = save_dir or os.getcwd()
    os.makedirs(save_dir, exist_ok=True)
    file_prefix = (
        f"{protein_file.file_prefix}_{ligand_file.file_prefix}"
        if ligand_file
        else protein_file.file_prefix
    )
    leap_inputfile = _create_leap_inputfile(
        protein_file_path=protein_file,
        ligand_file_path=ligand_file,
        frcmod_file_path=frcmod_file,
        file_prefix=file_prefix,
        box_size=box_size,
        box_type=box_type,
        solvatebox=solvatebox,
        save_dir=save_dir,
    )
    com_topfile_path = os.path.join(save_dir, f"{file_prefix}_com.prmtop")
    pro_topfile_path = os.path.join(save_dir, f"{file_prefix}_pro.prmtop")
    lig_topfile_path = os.path.join(save_dir, f"{file_prefix}_lig.prmtop")
    comsolvate_topfile_path = os.path.join(save_dir, f"{file_prefix}_comsolvate.prmtop")
    comsolvate_crdfile_path = os.path.join(save_dir, f"{file_prefix}_comsolvate.inpcrd")
    comsolvate_pdbfile_path = os.path.join(save_dir, f"{file_prefix}_comsolvate.pdb")
    with ChDir(save_dir):
        shell_run(f"tleap -f {leap_inputfile.file_path}")
        write_file(com_topfile_path, _get_strip_prmtop_str(File(com_topfile_path), DEFAULT_STRIP_MASKS))
        write_file(pro_topfile_path, _get_strip_prmtop_str(File(pro_topfile_path), DEFAULT_STRIP_MASKS))
        if ligand_file:
            write_file(lig_topfile_path, _get_strip_prmtop_str(File(lig_topfile_path), DEFAULT_STRIP_MASKS))

    shell_run(
        f"ambpdb -p {comsolvate_topfile_path} < {comsolvate_crdfile_path} > {comsolvate_pdbfile_path}"
    )
    return (
        File(comsolvate_topfile_path),
        File(comsolvate_crdfile_path),
        File(comsolvate_pdbfile_path),
        File(com_topfile_path),
        File(pro_topfile_path),
        File(lig_topfile_path) if ligand_file else None,
    )


class BaseProcess:
    """Base class for molecular dynamics simulation processes.

    Provides common functionality for parsing MD input files and managing
    simulation process parameters. This is an abstract base class that
    defines the interface for all MD simulation processes.

    Attributes:
        input_file (File): MD input file containing simulation parameters.
        process_name (str): Name identifier for this simulation process.
        cfg (list): Parsed configuration from the input file.
    """

    def __init__(self, input_file: File, process_name: str, **kwargs) -> None:
        """Initializes the base MD process.

        Args:
            input_file (File): Amber MD input file with simulation parameters.
            process_name (str): Unique name for this simulation process.
            **kwargs: Additional attributes to set on the process instance.

        Note:
            The input file is automatically parsed to extract configuration
            sections and parameters. Additional keyword arguments are set
            as instance attributes for process customization.
        """
        self.input_file = File(input_file)
        self.name = self.process_name = process_name
        self.cfg = self._get_input_config(input_file.file_path)
        for k, v in kwargs.items():
            setattr(self, k, v)

    @staticmethod
    def _get_input_config(input_file: str) -> list:
        """Parses Amber MD input file to extract configuration sections.

        Analyzes an Amber MD input file and extracts all configuration sections
        (namelist sections like &cntrl, &ewald, etc.) along with their parameters.

        Args:
            input_file (str): Path to the Amber MD input file to parse.

        Returns:
            list: List of dictionaries, each containing a configuration section with:
                - _index: Section order in the file
                - _type: Section type (e.g., 'cntrl', 'ewald')
                - parameter_name: parameter_value pairs from the section

        Note:
            The parser handles Amber's namelist format where sections start with
            & and end with /, and parameters are specified as key=value pairs.
            Empty sections and END markers are automatically skipped.
        """
        output_list = []
        type_pattern = r"&(.*?)\n"
        config_pattern = "(?<=&{_type}\n).*"
        item_pattern = r"(\w+)=(.*?)(?=(?:,\s*\w+=|$))"

        config_list = read_file(input_file).split("/")
        for index, config in enumerate(config_list):
            config = config.strip()
            if not config or config == "END":
                continue
            _type = re.findall(type_pattern, config)[0]
            config_dict = {"_index": index, "_type": _type}
            _config_list = re.findall(config_pattern.format(_type=_type), config, re.S)[0].replace(
                "\n", ""
            )
            for k, v in re.findall(item_pattern, _config_list):
                config_dict[k] = v
            output_list.append(config_dict)
        return output_list

    def run(self, **kwargs) -> None:
        """Executes the simulation process.

        Args:
            **kwargs: Process-specific arguments required for execution.

        Raises:
            NotImplementedError: This method must be implemented by subclasses.

        Note:
            Subclasses must override this method to implement the specific
            execution logic for their simulation type.
        """
        raise NotImplementedError


class MDProcess(BaseProcess):
    """Molecular dynamics simulation process implementation.

    Handles execution of Amber MD simulations including minimization, heating,
    equilibration, and production runs. Automatically parses input parameters
    and manages file paths for simulation output.

    Attributes:
        control_cfg (dict): Parsed &cntrl section from input file.
        is_minimize (bool): Whether this is a minimization process.
        is_restrained (bool): Whether restraints are applied.
        is_production (bool): Whether this is a production run.
        total_step (int): Total number of simulation steps.
        step_size (float): Time step size in picoseconds.
        cmd (str): Command string for simulation execution.
        topology_file (File): Topology file for the system.
        inpcrd_file (File): Input coordinate file.
        mdout_file_path (str): Path to simulation output file.
        mdcrd_file_path (str): Path to trajectory coordinate file.
        mdrst_file_path (str): Path to restart coordinate file.
    """

    def __init__(
        self,
        input_file: File,
        process_name: str,
        use_gpu: bool = True,
        cpu_num: int = None,
        **kwargs,
    ) -> None:
        """Initializes the MD process with input parameters.

        Args:
            input_file (File): Amber MD input file with simulation parameters.
            process_name (str): Unique name for this simulation process.
            use_gpu (bool, optional): Whether to use GPU acceleration. Defaults to True.
            cpu_num (int, optional): Number of CPU cores to use if not using GPU. Defaults to half of available CPUs.
            **kwargs: Additional attributes for process customization.

        Note:
            Automatically determines simulation type (minimization vs dynamics)
            based on the 'imin' parameter in the input file. Step size and
            total steps are extracted from appropriate input parameters.
        """
        super().__init__(input_file, process_name, **kwargs)
        self.control_cfg = [cfg for cfg in self.cfg if cfg["_type"] == "cntrl"][0]
        self.is_minimize = bool(int(self.control_cfg.get("imin") or 0))
        self.is_restrained = bool(int(self.control_cfg.get("ntr") or 0))
        self.is_production = False
        self.total_step = int(self.control_cfg.get("nstlim") or self.control_cfg.get("maxcyc"))
        self.step_size = float(self.control_cfg.get("dt")) if not self.is_minimize else 1
        self.use_gpu = use_gpu
        self.cpu_num = cpu_num or CPU_NUM
        self.cmd = ""
        self.topology_file = None
        self.inpcrd_file = None
        self.mdout_file_path = None
        self.mdcrd_file_path = None
        self.mdrst_file_path = None

    def run(
        self,
        topology_file: File,
        inpcrd_file: File,
        reference_file: File | str = None,
        save_dir: str = None,
        nohup: bool = False,
    ) -> "MDProcess":
        """Executes the molecular dynamics simulation.

        Runs the MD simulation using either PMEMD (for dynamics) or SANDER (for minimization)
        with the specified topology and coordinate files. Handles output file management
        and progress tracking.

        Args:
            topology_file (File): Amber topology file (.prmtop) for the system.
            inpcrd_file (File): Input coordinate file (.inpcrd or .rst).
            reference_file (File | str, optional): Reference structure for restrained simulations.
                Defaults to None.
            save_dir (str, optional): Directory to save simulation outputs.
                Defaults to current working directory.
            nohup (bool, optional): Whether to run in background with progress tracking.
                Defaults to False.

        Returns:
            MDProcess: Self-reference with updated file paths and execution status.

        Note:
            The method automatically:
            - Selects appropriate MD engine (PMEMD/SANDER) based on simulation type
            - Generates output file paths based on process name
            - Logs simulation parameters (steps, time step, total time)
            - Provides real-time progress tracking for background runs
            - Handles reference files for restrained simulations

        Raises:
            RuntimeError: If the simulation execution fails.
            FileNotFoundError: If required input files don't exist.
        """
        self.topology_file = topology_file
        self.inpcrd_file = inpcrd_file
        save_dir = save_dir or os.getcwd()
        os.makedirs(save_dir, exist_ok=True)

        self.mdout_file_path = os.path.join(save_dir, f"{self.process_name}.out")
        self.mdcrd_file_path = os.path.join(save_dir, f"{self.process_name}.nc")
        self.mdrst_file_path = os.path.join(save_dir, f"{self.process_name}.rst")
        simulation_time = self.total_step * self.step_size / 1000
        pmemd = PMEMD if self.use_gpu else SANDER.format(cpu_num=self.cpu_num)

        if not self.is_minimize:
            logger.info(f"Simulation total steps: {self.total_step}")
            logger.info(f"Simulation step size: {self.step_size} ps")
            logger.info(f"Simulation total time: {simulation_time} ns")

        self.cmd = f"{pmemd} -O"
        self.cmd += f" -i {self.input_file.file_path}"
        self.cmd += f" -o {self.mdout_file_path}"
        self.cmd += f" -c {self.inpcrd_file.file_path}"
        self.cmd += f" -p {self.topology_file.file_path}"
        self.cmd += f" -r {self.mdrst_file_path}"

        if not self.is_minimize:
            self.cmd += f" -x {self.mdcrd_file_path}"
        else:
            if reference_file is not None:
                self.cmd += f" -ref {reference_file.file_path}"
        if nohup:
            self.cmd = f"nohup {self.cmd} > /dev/null 2>&1 &"

        logger.debug(f"Process {self.process_name} running command: {self.cmd}")
        logger.info(f"Running process {self.process_name}...")
        shell_run(self.cmd)
        if nohup:
            _trace_progress(self.mdout_file_path, self.total_step)
        logger.info(f"Process {self.process_name} finished.")
        return self


class MinimizeProcess(MDProcess):
    """Energy minimization simulation process.

    Specialized MD process for performing energy minimization using SANDER.
    Inherits all functionality from MDProcess but is specifically configured
    for minimization calculations.

    Note:
        Minimization processes use SANDER engine and don't generate trajectory
        files. The process automatically detects minimization mode from the
        input file's 'imin=1' parameter.
    """

    def __init__(
        self,
        input_file: File,
        process_name: str,
        use_gpu: bool = False,
        cpu_num: int = None,
        **kwargs,
    ) -> None:
        """Initializes the minimization process.

        Args:
            input_file (File): Amber minimization input file with imin=1.
            process_name (str): Unique name for this minimization process.
            use_gpu (bool, optional): Whether to use GPU acceleration. Defaults to False.
            cpu_num (int, optional): Number of CPUs to use if not using GPU. Defaults to half of available CPUs.
            **kwargs: Additional attributes for process customization.
        """
        super().__init__(input_file, process_name, use_gpu=use_gpu, cpu_num=cpu_num, **kwargs)


class NPTProcess(MDProcess):
    """Isothermal-isobaric (NPT) molecular dynamics process.

    Performs MD simulation in the NPT ensemble where temperature and pressure
    are held constant. Commonly used for equilibration and production phases
    to maintain realistic density and pressure conditions.

    Attributes:
        is_production (bool): Whether this is a production run requiring
            background execution with progress tracking.
    """

    def __init__(
        self, input_file: File, process_name: str, is_production: bool = False, **kwargs
    ) -> None:
        """Initializes the NPT simulation process.

        Args:
            input_file (File): Amber MD input file configured for NPT simulation.
            process_name (str): Unique name for this NPT process.
            is_production (bool, optional): Whether this is a production run.
                Production runs use background execution. Defaults to False.
            **kwargs: Additional attributes for process customization.

        Note:
            NPT simulations require appropriate barostat and thermostat settings
            in the input file (ntp>0, ntt>0). Production runs automatically
            use nohup execution for long simulations.
        """
        super().__init__(input_file, process_name, **kwargs)
        self.is_production = is_production


class NVTProcess(MDProcess):
    """Canonical (NVT) molecular dynamics process.

    Performs MD simulation in the NVT ensemble where temperature and volume
    are held constant. Typically used for heating phases and some equilibration
    steps before switching to NPT conditions.

    Attributes:
        is_production (bool): Whether this is a production run requiring
            background execution with progress tracking.
    """

    def __init__(
        self, input_file: File, process_name: str, is_production: bool = False, **kwargs
    ) -> None:
        """Initializes the NVT simulation process.

        Args:
            input_file (File): Amber MD input file configured for NVT simulation.
            process_name (str): Unique name for this NVT process.
            is_production (bool, optional): Whether this is a production run.
                Production runs use background execution. Defaults to False.
            **kwargs: Additional attributes for process customization.

        Note:
            NVT simulations require thermostat settings in the input file (ntt>0)
            but no barostat (ntp=0). Commonly used for gradual heating from
            0 K to target temperature.
        """
        super().__init__(input_file, process_name, **kwargs)
        self.is_production = is_production


@timeit
def run_simulation(
    comsolvate_topfile: File,
    comsolvate_crdfile: File,
    process_list: list[MDProcess],
    save_dir: str = None,
) -> MDProcess:
    """Executes a series of molecular dynamics simulation processes.

    Runs a complete MD simulation workflow by sequentially executing a list
    of simulation processes (minimization, heating, equilibration, production).
    Each process uses the output coordinates from the previous step as input.

    Args:
        comsolvate_topfile (File): Amber topology file for the solvated system.
        comsolvate_crdfile (File): Initial coordinate file for the solvated system.
        process_list (list[MDProcess]): List of MD processes to execute in order.
            Each process contains its own input parameters and settings.
        save_dir (str, optional): Base directory to save all simulation outputs.
            Defaults to current working directory.

    Note:
        The function automatically:
        - Creates separate subdirectories for each simulation process
        - Chains coordinates between processes (output → input)
        - Handles restrained simulations with appropriate reference structures
        - Uses nohup for production runs to allow background execution
        - Provides progress tracking for long simulations

    Returns:
        MDProcess: The final completed simulation process in the workflow.

    Raises:
        RuntimeError: If any simulation process fails.
        FileNotFoundError: If input topology or coordinate files don't exist.
    """
    save_dir = save_dir or os.getcwd()
    os.makedirs(save_dir, exist_ok=True)

    for index, process in enumerate(process_list):
        process_output_dir = os.path.join(save_dir, process.process_name)
        os.makedirs(process_output_dir, exist_ok=True)
        inpcrd_file = comsolvate_crdfile if index == 0 else File(finished_process.mdrst_file_path)
        reference_file = inpcrd_file if process.is_restrained else None
        nohup = True if process.is_production else False
        finished_process = process.run(
            topology_file=comsolvate_topfile,
            inpcrd_file=inpcrd_file,
            save_dir=process_output_dir,
            reference_file=reference_file,
            nohup=nohup,
        )
    return finished_process


def create_energy_inputfile(
    job_type: str,
    startframe: int,
    endframe: int,
    step_size: int,
    decomp: bool = False,
    save_dir: str = None,
) -> File:
    """Creates MM-PBSA/MM-GBSA input file for binding energy calculations.

    Generates an input file for Amber's MM-PBSA or MM-GBSA calculations to
    analyze protein-ligand binding energies from MD trajectory data.

    Args:
        job_type (str): Type of energy calculation. Valid options:
            - "pb/gb": Both Poisson-Boltzmann and Generalized Born calculations
            - "gb": Generalized Born calculation only
            - "nmode": Normal mode entropy calculation
        startframe (int): Starting frame number from trajectory (1-indexed).
        endframe (int): Ending frame number from trajectory (1-indexed).
        step_size (int): Frame interval for analysis (e.g., 10 = every 10th frame).
        decomp (bool, optional): Whether to perform per-residue energy decomposition.
            Defaults to False.
        save_dir (str, optional): Directory to save the input file.
            Defaults to current working directory.

    Returns:
        File: Generated MM-PBSA/MM-GBSA input file.

    Raises:
        ValueError: If an invalid job_type is specified.

    Note:
        The function automatically configures appropriate calculation sections
        based on the job type. Energy decomposition provides detailed per-residue
        contributions but significantly increases computation time.
    """
    save_dir = save_dir or os.getcwd()
    os.makedirs(save_dir, exist_ok=True)
    mmgbsa_input = MMGBSAInput(start_frame=startframe, end_frame=endframe, step_size=step_size)
    mmgbsa_input.add_general()
    match job_type.lower():
        case "pb/gb":
            mmgbsa_input.add_pb()
            mmgbsa_input.add_gb()
            if decomp:
                mmgbsa_input.add_decomp()
        case "gb":
            mmgbsa_input.add_gb()
            if decomp:
                mmgbsa_input.add_decomp()
        case "nmode":
            mmgbsa_input.add_nmode()
        case _:
            raise ValueError(f"Invalid job_type: {job_type}")
    input_file_path = os.path.join(save_dir, f"{job_type.replace('/', '_')}.in")
    mmgbsa_input.save(input_file_path)
    return File(input_file_path)


@timeit
def run_energy_calculation(
    input_file: File,
    comsolvate_topfile: File,
    com_topfile: File,
    receptor_topfile: File,
    ligand_topfile: File,
    traj_file: File,
    save_dir: str = None,
    energy_decompose: bool = False,
    cpu_num: int = None,
) -> tuple[File, File]:
    """Executes MM-PBSA/MM-GBSA binding energy calculations.

    Performs comprehensive binding energy analysis using Amber's MMPBSA.py.MPI
    program with parallel processing. Calculates binding free energies and
    optionally performs per-residue energy decomposition.

    Args:
        input_file (File): MM-PBSA input file with calculation parameters.
        comsolvate_topfile (File): Topology file for the complete solvated system.
        com_topfile (File): Topology file for the protein-ligand complex (no solvent).
        receptor_topfile (File): Topology file for the receptor/protein only.
        ligand_topfile (File): Topology file for the ligand only.
        traj_file (File): MD trajectory file (.nc or .dcd format).
        save_dir (str, optional): Directory to save calculation results.
            Defaults to current working directory.
        energy_decompose (bool, optional): Whether per-residue decomposition was requested.
            Defaults to False.
        cpu_num (int, optional): Number of CPU cores for parallel calculation.
            Defaults to system CPU_NUM constant.

    Returns:
        tuple[File, File]: Tuple containing (energy_results, decomp_results):
            - energy_results: CSV file with binding energy results
            - decomp_results: CSV file with per-residue decomposition (if requested)

    Raises:
        RuntimeError: If MPI is not available or calculation fails.
        FileNotFoundError: If required input files don't exist.

    Note:
        The calculation requires:
        - MPI installation for parallel processing
        - Proper trajectory format compatibility with topology files
        - Sufficient memory for large trajectory analysis
        - All topology files must be consistent with the trajectory
    """
    save_dir = save_dir or os.getcwd()
    os.makedirs(save_dir, exist_ok=True)
    output_file = File(
        os.path.join(save_dir, f"{input_file.file_prefix}_FINAL_RESULTS_MMPBSA.csv"), exist=False
    )
    log_file = File(os.path.join(save_dir, f"{input_file.file_prefix}.log"), exist=False)
    decom_output_file = File(
        os.path.join(save_dir, f"{input_file.file_prefix}_FINAL_RESULTS_DECOMP.csv"), exist=False
    )
    logger.info("Checking MPI available...")
    if not is_mpirun_available():
        raise RuntimeError(
            "MPI is not available. Please install MPI/mpi4py to run energy calculations."
        )
    cpu_num = cpu_num or CPU_NUM
    energy_cmd = f"mpirun -np {cpu_num} --host localhost:{cpu_num} MMPBSA.py.MPI -O "
    energy_cmd += f"-i {input_file.file_path} "
    energy_cmd += f"-sp {comsolvate_topfile.file_path} "
    energy_cmd += f"-cp {com_topfile.file_path} "
    energy_cmd += f"-rp {receptor_topfile.file_path} "
    energy_cmd += f"-lp {ligand_topfile.file_path} "
    energy_cmd += f"-y {traj_file.file_path} "
    energy_cmd += f"-o {output_file.file_path} "
    if energy_decompose:
        energy_cmd += f"-do {decom_output_file.file_path} "
    energy_cmd += f"> {log_file.file_path} 2>&1"
    logger.info("Running Energy Calculation...")
    with ChDir():
        shell_run(energy_cmd)
    logger.info("Energy Calculation Finished.")
    return File(output_file), File(decom_output_file, exist=energy_decompose)
