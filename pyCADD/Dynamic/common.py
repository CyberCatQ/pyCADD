"""High-level molecular dynamics simulation workflow management.

This module provides the Processor class which serves as the main interface for
managing complete molecular dynamics simulation workflows including system preparation,
simulation execution, and trajectory analysis.
"""

import logging
import os
from typing import Literal

import pytraj as pt

from pyCADD.Dynamic import analysis, core
from pyCADD.Dynamic.constant import CPU_NUM
from pyCADD.Dynamic.core import MDProcess, MinimizeProcess, NPTProcess, NVTProcess
from pyCADD.Dynamic.template import (
    HeatInput,
    MinimizeInput,
    NPTInput,
    NVTInput,
    RestrainedMinimizeInput,
)
from pyCADD.Dynamic.utils import _get_water_resnum
from pyCADD.utils.common import File
from pyCADD.utils.tool import is_pmemd_cuda_available, makedirs_from_list

logger = logging.getLogger(__name__)


class Processor:
    """High-level molecular dynamics simulation workflow manager.

    The Processor class provides a complete interface for managing molecular dynamics
    simulations from system preparation through analysis.

    Attributes:
        work_dir (str): Main working directory for the simulation project.
        apo (bool): Whether this is an apo (protein-only) system.
        pro_file (File): Prepared protein structure file.
        mol2_file (File): Ligand MOL2 file (None for apo systems).
        frcmod_file (File): Force field parameter file for ligand.
        comsolvate_topfile (File): Solvated system topology file.
        comsolvate_crdfile (File): Solvated system coordinate file.
        comsolvate_pdbfile (File): Solvated system PDB file.
        com_topfile (File): Complex topology file.
        pro_topfile (File): Protein topology file.
        lig_topfile (File): Ligand topology file.
        md_process_list (list): List of MD simulation processes to execute.
    """

    def __init__(self, work_dir: str = None, apo: bool = False):
        """Initializes the MD simulation processor.

        Sets up the directory structure and initializes attributes for managing
        a molecular dynamics simulation workflow.

        Args:
            work_dir (str, optional): Main working directory for simulation files.
                Defaults to current working directory.
            apo (bool, optional): Whether this is an apo (protein-only) system.
                Defaults to False (protein-ligand system).

        Note:
            Automatically creates the following subdirectory structure:
            - protein/: For prepared protein files
            - molecule/: For ligand preparation files
            - leap/: For system building files
            - input_file/: For MD input files
            - md_result/: For simulation output files
            - md_analysis/: For trajectory analysis results
        """
        self.work_dir = work_dir or os.getcwd()
        self.pro_dir = os.path.join(self.work_dir, "protein")
        self.mol_dir = os.path.join(self.work_dir, "molecule")
        self.leap_dir = os.path.join(self.work_dir, "leap")
        self.input_file_dir = os.path.join(self.work_dir, "input_file")
        self.md_result_dir = os.path.join(self.work_dir, "md_result")
        self.analysis_result_dir = os.path.join(self.work_dir, "md_analysis")
        makedirs_from_list(self._required_dirs)

        self.pro_file = None
        self.mol2_file = None
        self.frcmod_file = None
        self.comsolvate_crdfile = None
        self.comsolvate_topfile = None
        self.comsolvate_pdbfile = None
        self.com_topfile = None
        self.pro_topfile = None
        self.lig_topfile = None
        if apo:
            logger.info("Initializing Processor for apo system.")
        self.apo = apo
        self.md_process_list = []

    @property
    def _required_dirs(self):
        return [
            self.pro_dir,
            self.mol_dir,
            self.leap_dir,
            self.input_file_dir,
            self.md_result_dir,
            self.analysis_result_dir,
        ]

    def protein_prepare(
        self, protein_file_path: str | File, keep_water: bool = False, overwrite: bool = False
    ) -> None:
        """Prepare protein structure for molecular dynamics simulation.

        Processes the input protein structure using pdb4amber to clean up the
        structure, add missing atoms, and prepare it for use with Amber force fields.

        Args:
            protein_file_path: Path to the input protein structure file (PDB format).
            keep_water: Whether to retain crystal water molecules in the structure.
                Defaults to False.
            overwrite: Whether to overwrite existing prepared files.
                Defaults to False.

        Note:
            The prepared protein file is automatically saved in the protein/
            subdirectory and stored in self.pro_file for subsequent use.
        """
        protein_file = File(protein_file_path)
        self.pro_file = core.protein_prepare(
            protein_file, save_dir=self.pro_dir, keep_water=keep_water, overwrite=overwrite
        )
        logger.info(f"Protein file has been saved: {self.pro_file.file_path}")

    def molecule_prepare(
        self,
        molecule_file_path: str | File,
        charge: int = 0,
        multiplicity: int = 1,
        dft: str = "B3LYP",
        basis_set: str = "6-31g*",
        cpu_num: int = 4,
        mem_use: str = "8GB",
        charge_method: Literal["resp", "bcc", "resp2"] = "resp",
        solvent: str = "water",
        delta: float = 0.5,
        overwrite: bool = False,
    ) -> None:
        """Prepare ligand molecule for molecular dynamics simulation.

        Performs complete ligand parameterization including geometry optimization,
        charge calculation, and force field parameter generation using antechamber
        and related Amber tools.

        Args:
            molecule_file_path: Path to the input ligand structure file.
            charge: Net charge of the ligand molecule. Defaults to 0.
            multiplicity: Spin multiplicity of the ligand. Defaults to 1.
            dft: DFT method for quantum chemistry calculations. Defaults to "B3LYP".
            basis_set: Basis set for quantum chemistry calculations. Defaults to "6-31g*".
            cpu_num: Number of CPU cores to use. Defaults to 4.
            mem_use: Memory allocation for quantum calculations. Defaults to "8GB".
            charge_method: Method for partial charge calculation ("resp", "bcc", "resp2").
                Defaults to "resp".
            solvent: Solvent model for charge calculation. Defaults to "water".
            delta: Grid spacing for RESP charge fitting. Defaults to 0.5.
            overwrite: Whether to overwrite existing prepared files. Defaults to False.

        Note:
            The prepared MOL2 and FRCMOD files are automatically saved in the
            molecule/ subdirectory and stored in self.mol2_file and self.frcmod_file.
        """

        cpu_num = cpu_num or CPU_NUM
        molecule_file = File(molecule_file_path)
        self.mol2_file, self.frcmod_file = core.molecule_prepare(
            ligand_file=molecule_file,
            charge_method=charge_method,
            charge=charge,
            multiplicity=multiplicity,
            dft=dft,
            basis_set=basis_set,
            cpu_num=cpu_num,
            mem_use=mem_use,
            solvent=solvent,
            delta=delta,
            overwrite=overwrite,
            save_dir=self.mol_dir,
        )
        logger.info(f"Molecule file has been saved: {self.mol2_file.file_path}")
        logger.info(f"Frcmod file has been saved: {self.frcmod_file.file_path}")

    def load_processed_profile(self, pro_file_path: str) -> None:
        """Load a previously prepared protein file.

        Args:
            pro_file_path: Path to the prepared protein structure file.
        """
        self.pro_file = File(pro_file_path)

    def load_processed_molfile(self, mol_file_path: str) -> None:
        """Load a previously prepared ligand MOL2 file.

        Args:
            mol_file_path: Path to the prepared ligand MOL2 file.
        """
        self.mol2_file = File(mol_file_path)

    def load_frcmod_file(self, frcmod_file_path: str) -> None:
        """Load a ligand force field parameter file.

        Args:
            frcmod_file_path: Path to the ligand FRCMOD parameter file.
        """
        self.frcmod_file = File(frcmod_file_path)

    def leap_prepare(
        self,
        box_size: float = 10.0,
        box_type: str = "TIP3PBOX",
        solvatebox: Literal["solvatebox", "solvateOct"] = "solvatebox",
    ) -> None:
        """Prepare solvated molecular system using tLEaP.

        Builds the complete simulation system by combining protein and ligand
        (if applicable), adding solvent and counter-ions using Amber's tLEaP program.

        Args:
            box_size: Distance from solute to box edge in Angstroms. Defaults to 10.0.
            box_type: Type of water box ("TIP3PBOX", "TIP4PBOX", etc.). Defaults to "TIP3PBOX".
            solvatebox: Solvation command type ("solvatebox" or "solvateOct").
                Defaults to "solvatebox".

        Raises:
            RuntimeError: If required protein or ligand files haven't been prepared.

        Note:
            Automatically generates topology (.top), coordinate (.crd), and PDB files
            for the solvated system in the leap/ subdirectory.
        """
        box_type = box_type.upper()
        _check_dirs = (
            [self.mol2_file, self.frcmod_file, self.pro_file] if not self.apo else [self.pro_file]
        )
        if not all(_check_dirs):
            raise RuntimeError("Preparing protein or molecules before LEaP has not been done.")
        logger.info("Preparing LEaP files...")
        logger.info(f"Creating {box_type} water box with '{solvatebox}'")
        logger.info(f"{box_type} water box size: {box_size} Angstroms")

        (
            self.comsolvate_topfile,
            self.comsolvate_crdfile,
            self.comsolvate_pdbfile,
            self.com_topfile,
            self.pro_topfile,
            self.lig_topfile,
        ) = core.leap_prepare(
            ligand_file=self.mol2_file,
            frcmod_file=self.frcmod_file,
            protein_file=self.pro_file,
            box_size=box_size,
            box_type=box_type,
            solvatebox=solvatebox,
            save_dir=self.leap_dir,
        )

        logger.info(f"LEaP files have been saved in {self.leap_dir}")

    def set_prepared_file(self, file_path: str, file_type: str) -> None:
        """Set a pre-prepared file for use in the simulation workflow.

        Allows manual specification of prepared files, useful for resuming workflows
        or using externally prepared files.

        Args:
            file_path: Path to the prepared file.
            file_type: Type of file to set. Valid options:
                - "pro": Prepared protein file
                - "mol2": Prepared ligand MOL2 file
                - "frcmod": Ligand force field parameter file
                - "com_pdb": Solvated complex PDB file
                - "com_top": Solvated complex topology file
                - "com_crd": Solvated complex coordinate file

        Raises:
            RuntimeError: If an invalid file_type is specified.
        """
        _file = File(file_path)
        match file_type:
            case "pro":
                logger.info(f"Setting prepared protein file: {file_path}")
                self.pro_file = _file
            case "mol2":
                logger.info(f"Setting prepared mol2 file: {file_path}")
                self.mol2_file = _file
            case "frcmod":
                logger.info(f"Setting prepared frcmod file: {file_path}")
                self.frcmod_file = _file
            case "com_pdb":
                logger.info(f"Setting prepared solvated complex pdb file: {file_path}")
                self.comsolvate_pdbfile = _file
            case "com_top":
                logger.info(f"Setting prepared solvated complex topology file: {file_path}")
                self.comsolvate_topfile = _file
            case "com_crd":
                logger.info(f"Setting prepared solvated complex coordinate file: {file_path}")
                self.comsolvate_crdfile = _file
            case _:
                raise RuntimeError(f"{file_type} is not a valid file type.")

    def get_water_resnum(self) -> list:
        """Get list of water residue numbers from the solvated system.

        Returns:
            list: List of residue numbers corresponding to water molecules.

        Raises:
            ValueError: If the solvated complex PDB file hasn't been prepared or loaded.

        Note:
            Useful for setting up restraints or analysis that excludes water molecules.
        """
        if self.comsolvate_pdbfile is None:
            raise ValueError("Solvated complex pdb file has not been prepared/load.")
        return _get_water_resnum(self.comsolvate_pdbfile)

    def create_minimize_input(
        self,
        maxcyc: int = 10000,
        ncyc: int = 5000,
        cut: float = 8.0,
        restraint: bool = False,
        restraint_mask: str = None,
        restraint_wt: float = 2.0,
        file_name: str = None,
        **kwargs,
    ) -> File:
        """Create energy minimization input file.

        Generates Amber input files for energy minimization using steepest descent
        followed by conjugate gradient algorithms, with optional positional restraints.

        Args:
            maxcyc: Maximum number of minimization cycles. Defaults to 10000.
            ncyc: Number of steepest descent cycles before switching to conjugate gradient.
                Defaults to 5000.
            cut: Nonbonded interaction cutoff distance in Angstroms. Defaults to 8.0.
            restraint: Whether to apply positional restraints. Defaults to False.
            restraint_mask: Amber mask syntax for restrained atoms. Required if restraint=True.
            restraint_wt: Restraint force constant in kcal/mol/Å². Defaults to 2.0.
            file_name: Output filename for the input file. Defaults to "minimize.in".
            **kwargs: Additional minimization parameters.

        Returns:
            File: The created minimization input file.

        Raises:
            RuntimeError: If restraint=True but no restraint_mask is provided.
        """
        if not restraint:
            constructor = MinimizeInput(maxcyc=maxcyc, ncyc=ncyc, cut=cut, **kwargs)
        else:
            if restraint_mask is None:
                raise RuntimeError("restraint_mask is required when restraint is True.")
            constructor = RestrainedMinimizeInput(
                maxcyc=maxcyc,
                ncyc=ncyc,
                cut=cut,
                restraintmask=restraint_mask,
                restraint_wt=restraint_wt,
                **kwargs,
            )
        file_name = file_name or "minimize.in"
        file_path = os.path.join(self.input_file_dir, file_name)
        constructor.save(file_path)
        logger.info(f'Minimize process input file created for {file_name.split(".")[0]}.')
        logger.info(f"{constructor.get_state_dict()}")
        return File(file_path)

    def create_heat_input(
        self,
        target_temp: float = 300.0,
        heat_step: int = 9000,
        total_step: int = 10000,
        step_length: float = 0.002,
        restraint_mask: str = None,
        restraint_wt: float = 2.0,
        file_name: str = None,
    ) -> File:
        """Create system heating input file.

        Generates Amber input file for gradual temperature increase from 0 K to
        target temperature using a controlled heating protocol with optional restraints.

        Args:
            target_temp: Target temperature in Kelvin. Defaults to 300.0.
            heat_step: Number of steps for temperature ramping phase. Defaults to 9000.
            total_step: Total simulation steps for heating protocol. Defaults to 10000.
            step_length: Time step length in picoseconds. Defaults to 0.002.
            restraint_mask: Optional Amber mask for restrained atoms.
            restraint_wt: Restraint force constant in kcal/mol/Å². Defaults to 2.0.
            file_name: Output filename for input file. Defaults to "heat.in".

        Returns:
            File: The created heating input file.

        Note:
            Uses a two-phase protocol: linear temperature ramp followed by
            temperature maintenance at the target value.
        """
        constructor = HeatInput(
            target_temp=target_temp,
            heat_step=heat_step,
            total_step=total_step,
            step_length=step_length,
            restraint_wt=restraint_wt,
            restraintmask=restraint_mask,
        )
        file_name = file_name or "heat.in"
        file_path = os.path.join(self.input_file_dir, file_name)
        constructor.save(file_path)
        logger.info(f'Heat process input file created for {file_name.split(".")[0]}')
        logger.info(f"{constructor.get_state_dict()}")
        return File(file_path)

    def create_nvt_input(
        self,
        init_temp: float = 300.0,
        total_step: int = 500000,
        step_length: float = 0.002,
        irest: int = 1,
        ntx: int = 5,
        file_name: str = None,
        **kwargs,
    ) -> File:
        """Create NVT (constant volume, temperature) simulation input file.

        Generates Amber input file for NVT ensemble molecular dynamics simulation
        with Langevin thermostat for temperature control.

        Args:
            init_temp: Simulation temperature in Kelvin. Defaults to 300.0.
            total_step: Total number of MD steps. Defaults to 500000 (1 ns at 2 fs steps).
            step_length: Time step length in picoseconds. Defaults to 0.002.
            irest: Restart flag (0=new simulation, 1=restart). Defaults to 1.
            ntx: Coordinate/velocity reading option. Defaults to 5.
            file_name: Output filename for input file. Defaults to "nvt.in".
            **kwargs: Additional NVT simulation parameters.

        Returns:
            File: The created NVT input file.

        Note:
            Typically used for equilibration after heating or as a production
            simulation for constant volume systems.
        """
        constructor = NVTInput(
            temp0=init_temp, nstlim=total_step, dt=step_length, irest=irest, ntx=ntx, **kwargs
        )
        file_name = file_name or "nvt.in"
        file_path = os.path.join(self.input_file_dir, file_name)
        constructor.save(file_path)
        logger.info(f'NVT process input file created for {file_name.split(".")[0]}')
        logger.info(f"{constructor.get_state_dict()}")
        return File(file_path)

    def create_npt_input(
        self,
        init_temp: float = 300.0,
        total_step: int = 50000000,
        step_length: float = 0.002,
        irest: int = 1,
        ntx: int = 5,
        taup: float = 2.0,
        file_name: str = None,
        **kwargs,
    ) -> File:
        """Create NPT (constant pressure, temperature) simulation input file.

        Generates Amber input file for NPT ensemble molecular dynamics simulation
        with both temperature and pressure control for production simulations.

        Args:
            init_temp: Simulation temperature in Kelvin. Defaults to 300.0.
            total_step: Total number of MD steps. Defaults to 50000000 (100 ns at 2 fs steps).
            step_length: Time step length in picoseconds. Defaults to 0.002.
            irest: Restart flag (0=new simulation, 1=restart). Defaults to 1.
            ntx: Coordinate/velocity reading option. Defaults to 5.
            taup: Pressure relaxation time in picoseconds. Defaults to 2.0.
            file_name: Output filename for input file. Defaults to "npt.in".
            **kwargs: Additional NPT simulation parameters.

        Returns:
            File: The created NPT input file.

        Note:
            Typically used for production simulations at physiological conditions
            (1 atm pressure, 300 K temperature).
        """
        constructor = NPTInput(
            temp0=init_temp,
            nstlim=total_step,
            dt=step_length,
            irest=irest,
            ntx=ntx,
            taup=taup,
            **kwargs,
        )
        file_name = file_name or "npt.in"
        file_path = os.path.join(self.input_file_dir, file_name)
        constructor.save(file_path)
        logger.info(f'NPT process input file created for {file_name.split(".in")[0]}')
        logger.info(f"{constructor.get_state_dict()}")
        return File(file_path)

    def add_process(
        self,
        input_file: str | File,
        process_name: str = None,
        _type: str = None,
        process: MDProcess = None,
        **kwargs,
    ) -> None:
        """Add a molecular dynamics process to the simulation workflow.

        Adds an MD simulation process (minimization, heating, equilibration, or production)
        to the execution queue with appropriate process type and parameters.

        Args:
            input_file: Path to the Amber input file for this process.
            process_name: Descriptive name for this process step.
            _type: Type of MD process ("minimize", "nvt", "npt"). Defaults to generic MDProcess.
            process: Custom MDProcess subclass to use instead of auto-detected type.
            **kwargs: Additional parameters passed to the process constructor.

        Note:
            Process type determines the appropriate MDProcess subclass:
            - "minimize": MinimizeProcess
            - "nvt": NVTProcess
            - "npt": NPTProcess
            - default: generic MDProcess
        """
        input_file = File(input_file)
        if _type == "minimize":
            process_class = MinimizeProcess
        elif _type == "nvt":
            process_class = NVTProcess
        elif _type == "npt":
            process_class = NPTProcess
        else:
            process_class = MDProcess
        if process is not None:
            process_class = process
        self.md_process_list.append(process_class(input_file, process_name, **kwargs))
        logger.info(f"Process {process_name} has been added to Workflow.")

    def add_minimize_process(
        self,
        maxcyc: int = 10000,
        ncyc: int = 5000,
        process_name: str = "minimize",
        restraint: bool = False,
        restraint_mask: str = None,
        use_gpu: bool = False,
        cpu_num: int = None,
        **kwargs,
    ) -> None:
        """Add energy minimization process to workflow.

        Creates a minimization input file and adds the corresponding process
        to the simulation workflow queue.

        Args:
            maxcyc: Maximum number of minimization cycles. Defaults to 10000.
            ncyc: Number of steepest descent cycles. Defaults to 5000.
            process_name: Name for this minimization step. Defaults to "minimize".
            restraint: Whether to apply positional restraints. Defaults to False.
            restraint_mask: Amber mask for restrained atoms if restraint=True.
            use_gpu: Whether to use GPU acceleration if available. Defaults to False.
            cpu_num: Number of CPU cores to use if not using GPU. Defaults to half of available CPUs.
            **kwargs: Additional minimization parameters.

        Note:
            Automatically creates input file and adds MinimizeProcess to workflow.
        """
        file_name = f"{process_name}.in"
        self.add_process(
            self.create_minimize_input(
                maxcyc=maxcyc,
                ncyc=ncyc,
                restraint=restraint,
                restraint_mask=restraint_mask,
                file_name=file_name,
                **kwargs,
            ),
            process_name,
            _type="minimize",
            use_gpu=use_gpu,
            cpu_num=cpu_num,
            **kwargs,
        )

    def add_nvt_process(
        self,
        total_step: int = 500000,
        step_length: float = 0.002,
        process_name: str = "nvt",
        is_production: bool = False,
        **kwargs,
    ):
        """Add NVT simulation process to workflow.

        Creates an NVT input file and adds the corresponding process
        to the simulation workflow queue.

        Args:
            total_step: Total number of MD steps. Defaults to 500000 (1 ns).
            step_length: Time step length in picoseconds. Defaults to 0.002.
            process_name: Name for this NVT step. Defaults to "nvt".
            is_production: Whether this is a production simulation. Defaults to False.
            **kwargs: Additional NVT simulation parameters.

        Note:
            Automatically creates input file and adds NVTProcess to workflow.
        """
        file_name = f"{process_name}.in"
        self.add_process(
            self.create_nvt_input(
                total_step=total_step, step_length=step_length, file_name=file_name, **kwargs
            ),
            process_name,
            _type="nvt",
            is_production=is_production,
            **kwargs,
        )

    def add_heat_process(
        self,
        target_temp: float = 300.0,
        heat_step: int = 9000,
        total_step: int = 10000,
        step_length: float = 0.002,
        process_name: str = "heat",
        **kwargs,
    ):
        """Add system heating process to workflow.

        Creates a heating input file and adds the corresponding process
        to the simulation workflow queue.

        Args:
            target_temp: Target temperature in Kelvin. Defaults to 300.0.
            heat_step: Number of steps for temperature ramping. Defaults to 9000.
            total_step: Total simulation steps. Defaults to 10000.
            step_length: Time step length in picoseconds. Defaults to 0.002.
            process_name: Name for this heating step. Defaults to "heat".
            **kwargs: Additional heating parameters.

        Note:
            Automatically creates input file and adds MDProcess to workflow.
        """
        file_name = f"{process_name}.in"
        self.add_process(
            self.create_heat_input(
                target_temp=target_temp,
                heat_step=heat_step,
                total_step=total_step,
                step_length=step_length,
                file_name=file_name,
            ),
            process_name,
            **kwargs,
        )

    def add_npt_process(
        self,
        total_step: int = 50000000,
        step_length: float = 0.002,
        process_name: str = "npt",
        is_production: bool = False,
        **kwargs,
    ):
        """Add NPT simulation process to workflow.

        Creates an NPT input file and adds the corresponding process
        to the simulation workflow queue.

        Args:
            total_step: Total number of MD steps. Defaults to 50000000 (100 ns).
            step_length: Time step length in picoseconds. Defaults to 0.002.
            process_name: Name for this NPT step. Defaults to "npt".
            is_production: Whether this is a production simulation. Defaults to False.
            **kwargs: Additional NPT simulation parameters.

        Note:
            Automatically creates input file and adds NPTProcess to workflow.
        """
        file_name = f"{process_name}.in"
        self.add_process(
            self.create_npt_input(
                total_step=total_step, step_length=step_length, file_name=file_name, **kwargs
            ),
            process_name,
            _type="npt",
            is_production=is_production,
            **kwargs,
        )


class Simulator:
    """GPU-accelerated molecular dynamics simulation executor.

    Manages the execution of MD simulation workflows on GPU hardware using
    AMBER's pmemd.cuda engine. Handles GPU device selection and process
    execution coordination.

    Attributes:
        processor: Associated Processor instance containing the workflow.
        comsolvate_topfile: System topology file for simulation.
        comsolvate_crdfile: System coordinate file for simulation.
        md_process_list: List of MD processes to execute.
        cuda_device: GPU device ID to use for calculations.
    """

    def __init__(self, processor: Processor) -> None:
        """Initialize Simulator with a configured Processor.

        Args:
            processor: Processor instance with prepared system files and
                process workflow.

        Note:
            Automatically inherits system files and process list from
            the provided processor instance.
        """
        self.processor = processor
        self.comsolvate_topfile = processor.comsolvate_topfile
        self.comsolvate_crdfile = processor.comsolvate_crdfile
        self.com_topfile = processor.com_topfile
        self.pro_topfile = processor.pro_topfile
        self.lig_topfile = processor.lig_topfile
        self.traj_file = None
        self.mdout_file = None
        self.md_process_list = processor.md_process_list
        self.cuda_device = 0

    def show_cuda_device(self) -> None:
        """Display available GPU devices using nvidia-smi.

        Prints GPU information including memory usage, temperature,
        and process information for device selection.
        """
        gpu_info = os.popen("nvidia-smi").read()
        print(gpu_info)

    def set_cuda_device(self, cuda_device: int) -> None:
        """Set the GPU device for CUDA calculations.

        Args:
            cuda_device: GPU device ID to use (0, 1, 2, etc.).

        Note:
            Sets the CUDA_VISIBLE_DEVICES environment variable to
            restrict CUDA operations to the specified device.
        """
        self.cuda_device = cuda_device
        os.environ["CUDA_VISIBLE_DEVICES"] = str(cuda_device)
        logger.info(f"Using GPU device {cuda_device}")

    def run_simulation(self, cuda_device: int = None) -> None:
        """Execute the complete MD simulation workflow on GPU.

        Runs all MD processes in the workflow using AMBER's pmemd.cuda
        engine with automatic file management and error checking.

        Args:
            cuda_device: Optional GPU device ID to use. If provided,
                overrides the current device setting.

        Raises:
            RuntimeError: If no MD processes have been added to workflow
                or if pmemd.cuda is not available.

        Note:
            Results are automatically saved in the md_result/ subdirectory
            with appropriate naming conventions for each process step.
        """
        if len(self.md_process_list) == 0:
            raise RuntimeError("No MD process has been added.")
        if not is_pmemd_cuda_available():
            raise RuntimeError("pmemd.cuda is not installed or not in PATH.")
        if cuda_device is not None:
            self.set_cuda_device(cuda_device)

        production_process = core.run_simulation(
            comsolvate_topfile=self.comsolvate_topfile,
            comsolvate_crdfile=self.comsolvate_crdfile,
            process_list=self.md_process_list,
            save_dir=self.processor.md_result_dir,
        )
        self.traj_file = File(production_process.mdcrd_file_path)
        self.mdout_file = File(production_process.mdout_file_path)
        logger.info(f"Simulation workflow finished.")


class Analyzer:
    """Molecular dynamics trajectory analysis interface.

    Provides comprehensive analysis capabilities for MD simulation trajectories
    including RMSD, RMSF, hydrogen bond analysis, and energy extraction using
    pytraj and custom analysis functions.

    Attributes:
        traj_file_path: Path to the MD trajectory file.
        traj_file: File object for trajectory.
        traj: Loaded trajectory object for analysis.
        top_file_path: Path to the topology file.
        top_file: File object for topology.
    """

    def __init__(
        self,
        traj_file: str | File,
        comsolvated_topfile: str | File,
        com_topfile: str | File = None,
        receptor_topfile: str | File = None,
        ligand_topfile: str | File = None,
        mdout_file: str | File = None,
        save_dir: str = None,
    ) -> None:
        """Initialize trajectory analyzer with file paths.

        Args:
            traj_file (str | File): Path to MD trajectory file (.nc, .mdcrd, etc.).
            comsolvated_topfile (str | File): Path to solvated system topology file.
            com_topfile (str | File): Path to complex (protein-ligand) topology file.
            receptor_topfile (str | File): Path to receptor-only topology file.
            ligand_topfile (str | File): Path to ligand-only topology file.
            mdout_file (str | File): Path to MD output file for energy extraction.
            save_dir: Directory for saving analysis results.

        Note:
            Multiple topology files support different analysis types
            (complex, apo, ligand-only). Trajectory and topology files
            are loaded on-demand when analysis methods are called.
        """
        self.traj_file = File(traj_file)
        self.top_file = File(comsolvated_topfile)
        self.load_traj(traj_file, comsolvated_topfile)

        self.mdout_file_path = File(mdout_file) if mdout_file is not None else None
        self.com_topfile = File(com_topfile) if com_topfile is not None else None
        self.recep_topfile = File(receptor_topfile) if receptor_topfile is not None else None
        self.apo = False if ligand_topfile is not None else True
        self.ligand_topfile = File(ligand_topfile) if ligand_topfile is not None else None
        self.save_dir = save_dir or os.getcwd()
        os.makedirs(self.save_dir, exist_ok=True)

        self.load_topfile(
            com_topfile=self.com_topfile,
            receptor_topfile=self.recep_topfile,
            ligand_topfile=self.ligand_topfile,
        )

        self.rmsd = None
        self.rmsf = None
        self.hbond = None
        self.distance = None
        self.angle = None

    def load_traj(self, traj_file: str | File, top_file: str | File) -> None:
        """Load trajectory and topology files for analysis.

        Args:
            traj_file (str | File): Path to the MD trajectory file.
            top_file (str | File): Path to the corresponding topology file.

        Note:
            Uses pytraj to load trajectory for efficient analysis operations.
            Prints trajectory information including frame count and atom count.
        """
        self.traj_file = File(traj_file)
        self.top_file = File(top_file)
        self.traj = pt.iterload(self.traj_file.file_path, self.top_file.file_path)
        logger.info(f"Trajectory file has been loaded: {self.traj_file.file_path}")
        logger.info(f"Topology file has been loaded: {self.top_file.file_path}")
        logger.info(f"Trajectory Info: {self.traj}")

    def load_mdout(self, mdout_file: str | File) -> None:
        """Load MD output file for energy analysis.

        Args:
            mdout_file (str | File): Path to the Amber MD output file containing
                energy information from the simulation.
        """
        self.mdout_file = File(mdout_file)
        logger.info(f"MD output file has been loaded: {self.mdout_file.file_path}")

    def load_topfile(
        self,
        comsolvated_topfile: str | File = None,
        com_topfile: str | File = None,
        receptor_topfile: str | File = None,
        ligand_topfile: str | File = None,
    ) -> None:
        """Load topology files for different system components.

        Args:
            comsolvated_topfile: Path to solvated complex topology.
            com_topfile: Path to complex (protein-ligand) topology.
            receptor_topfile: Path to receptor-only topology.
            ligand_topfile: Path to ligand-only topology.

        Note:
            Multiple topology files enable component-specific analysis
            (e.g., protein-only RMSD, ligand-only RMSF). Sets apo flag
            based on ligand topology availability.
        """
        if comsolvated_topfile is not None:
            self.top_file = File(comsolvated_topfile)
            logger.info(
                f"Solvated complex topology file has been loaded: {self.top_file.file_path}"
            )
        if com_topfile is not None:
            self.com_topfile = File(com_topfile)
            logger.info(f"Complex topology file has been loaded: {self.com_topfile.file_path}")
        if receptor_topfile is not None:
            self.recep_topfile = File(receptor_topfile)
            logger.info(f"Receptor topology file has been loaded: {self.recep_topfile.file_path}")
        if ligand_topfile is not None:
            self.ligand_topfile = File(ligand_topfile)
            logger.info(f"Ligand topology file has been loaded: {self.ligand_topfile.file_path}")
        else:
            logger.info("Ligand topology file is not provided. Perform analysis for Apo system.")

    def calc_rmsd(self, mask: str = "@CA", ref: int = 0, **kwargs) -> None:
        """Calculate Root Mean Square Deviation (RMSD) over trajectory.

        Computes RMSD values for selected atoms relative to a reference frame,
        providing insight into structural stability and conformational changes.

        Args:
            mask: Amber mask syntax for atom selection. Defaults to "@CA" (alpha carbons).
            ref: Reference frame number for RMSD calculation. Defaults to 0 (first frame).
            **kwargs: Additional parameters passed to analysis function.

        Note:
            Results are saved in rmsd/ subdirectory as CSV files and plots.
            Common masks: "@CA" (backbone), ":1-100" (residues), "!:WAT" (no water).
        """
        save_dir = os.path.join(self.save_dir, "rmsd")
        os.makedirs(save_dir, exist_ok=True)
        logger.info(f"Calculating RMSD...")
        logger.info(f"Amber mask: {mask}")
        logger.info(f"Reference frame: {ref}")
        self.rmsd = analysis._calc_rmsd(
            self.traj, mask=mask, reference=ref, save_dir=save_dir, **kwargs
        )

    def calc_rmsf(self, mask: str = "@CA", options: str = "byres", **kwargs) -> None:
        """Calculate Root Mean Square Fluctuation (RMSF) per residue/atom.

        Computes RMSF values to identify flexible regions and binding site
        dynamics in the protein structure.

        Args:
            mask: Amber mask syntax for atom selection. Defaults to "@CA".
            options: Analysis options ("byres" for per-residue, "byatom" for per-atom).
                Defaults to "byres".
            **kwargs: Additional parameters passed to analysis function.

        Note:
            Results are saved in rmsf/ subdirectory. "byres" option averages
            fluctuations over all atoms in each residue for cleaner plots.
        """
        save_dir = os.path.join(self.save_dir, "rmsf")
        os.makedirs(save_dir, exist_ok=True)
        logger.info(f"Calculating RMSF...")
        logger.info(f"Amber mask: {mask}")
        logger.info(f"Options: {options}")
        self.rmsf = analysis._calc_rmsf(
            self.traj, mask=mask, options=options, save_dir=save_dir, **kwargs
        )

    def calc_hbond(
        self,
        mask: str = ":*",
        distance: float = 3.0,
        angle: float = 135.0,
        options: str = None,
        **kwargs,
    ) -> None:
        """Calculate hydrogen bond analysis throughout trajectory.

        Identifies and tracks hydrogen bonds based on geometric criteria,
        providing insights into protein-ligand and internal protein interactions.

        Args:
            mask: Amber mask syntax for donor/acceptor selection. Defaults to ":*" (all).
            distance: Maximum donor-acceptor distance in Angstroms. Defaults to 3.0.
            angle: Minimum donor-H-acceptor angle in degrees. Defaults to 135.0.
            options: Additional options for hydrogen bond analysis.
            **kwargs: Additional parameters passed to analysis function.

        Note:
            Results include hydrogen bond occupancy, distances, and angles
            saved in hbond/ subdirectory with detailed statistics.
        """
        save_dir = os.path.join(self.save_dir, "hbond")
        os.makedirs(save_dir, exist_ok=True)
        logger.info(f"Calculating and tracing H-bonds...")
        self.hbond = analysis._calc_hbond(
            self.traj,
            mask=mask,
            distance=distance,
            angle=angle,
            options=options,
            save_dir=save_dir,
            **kwargs,
        )

    def trace_distance(self, mask: str, **kwargs) -> None:
        """Trace distance between two atom groups over trajectory.

        Monitors the distance between specified atoms/groups to analyze
        binding interactions, conformational changes, or domain movements.

        Args:
            mask: Amber mask syntax specifying two atom groups (e.g., ":1@CA :100@CA").
            **kwargs: Additional parameters passed to analysis function.

        Note:
            Results saved in distance/ subdirectory with time series plots
            and statistical analysis of distance fluctuations.
        """
        save_dir = os.path.join(self.save_dir, "distance")
        os.makedirs(save_dir, exist_ok=True)
        self.distance = analysis._trace_distance(
            self.traj, mask=mask, save=True, save_dir=save_dir, **kwargs
        )

    def trace_angle(self, mask: str, **kwargs) -> None:
        """Trace angle between three atom groups over trajectory.

        Monitors angular measurements to analyze conformational changes,
        binding geometry, or secondary structure dynamics.

        Args:
            mask: Amber mask syntax specifying three atom groups for angle calculation
                (e.g., ":1@CA :2@CA :3@CA").
            **kwargs: Additional parameters passed to analysis function.

        Note:
            Results saved in angle/ subdirectory with time series plots
            showing angle evolution and distribution analysis.
        """
        save_dir = os.path.join(self.save_dir, "angle")
        os.makedirs(save_dir, exist_ok=True)
        self.angle = analysis._trace_angle(
            self.traj, mask=mask, save=True, save_dir=save_dir, **kwargs
        )

    def extract_frame(self, frame: int, **kwargs) -> None:
        """Extract a single frame structure from trajectory.

        Saves the specified frame as a separate structure file for
        detailed analysis or visualization.

        Args:
            frame: Frame number to extract (0-indexed).
            **kwargs: Additional parameters passed to extraction function.

        Note:
            Frame structures saved in frame_structures/ subdirectory
            in PDB format for further analysis or visualization.
        """
        save_dir = os.path.join(self.save_dir, "frame_structures")
        os.makedirs(save_dir, exist_ok=True)
        self.frame = analysis._extract_frame(
            traj=self.traj, frame_indices=[frame], save_dir=save_dir, **kwargs
        )

    def extract_frames(self, start: int, end: int, **kwargs) -> None:
        """Extract a range of frame structures from trajectory.

        Saves multiple frames as separate structure files for
        ensemble analysis or representative structure selection.

        Args:
            start: Starting frame number (0-indexed, inclusive).
            end: Ending frame number (0-indexed, inclusive).
            **kwargs: Additional parameters passed to extraction function.

        Note:
            All frames in the specified range are saved in frame_structures/
            subdirectory with sequential numbering.
        """
        save_dir = os.path.join(self.save_dir, "frame_structures")
        os.makedirs(save_dir, exist_ok=True)
        indices = range(start, end + 1)
        self.frames = analysis._extract_frame(
            traj=self.traj, frame_indices=indices, save_dir=save_dir, **kwargs
        )

    def extract_lowest_energy_st(self) -> None:
        """Extract structure frame with lowest total energy.

        Analyzes the MD output file to identify the frame with minimum total energy
        and extracts it as a representative low-energy structure.

        Raises:
            ValueError: If MD output file hasn't been loaded.

        Note:
            Results include energy information and extracted structure saved
            in lowest_energy_structure/ subdirectory. Useful for identifying
            stable conformations or starting structures for further analysis.
        """
        save_dir = os.path.join(self.save_dir, "lowest_energy_structure")
        os.makedirs(save_dir, exist_ok=True)
        if self.mdout_file is None:
            raise ValueError("Mdout file (.out) is required for extracting lowest energy structure.")

        logger.info(f"Detecting lowest energy frame...")
        self.LE_frame, self.LE_time, self.LE_energy = analysis._get_lowest_energy_info(
            mdout_file=self.mdout_file, save_dir=save_dir
        )
        logger.info(f"Lowest energy frame: {self.LE_frame}")
        logger.info(f"Lowest energy time: {self.LE_time}")
        logger.info(f"Lowest energy: {self.LE_energy}")
        logger.info(f"Extracting lowest energy structure: Frame {self.LE_frame}...")
        analysis._extract_frame(self.traj, frame_indices=[self.LE_frame], save_dir=save_dir)

    def create_energy_inputfile(
        self,
        start_frame: int,
        end_frame: int,
        job_type: Literal["free", "entropy", "decomp"],
        method: Literal["pb/gbsa", "gbsa"] = None,
        step_size: int = 10,
    ) -> None:
        """Create input file for MM/PBSA or MM/GBSA energy calculations.

        Generates input files for free energy analysis, entropy calculations,
        or energy decomposition using AMBER's MMPBSA.py tool.

        Args:
            start_frame: First frame for energy analysis (1-indexed).
            end_frame: Last frame for energy analysis (inclusive).
            job_type: Type of energy calculation:
                - "free": Free energy calculation (ΔG)
                - "entropy": Entropy contribution calculation
                - "decomp": Per-residue energy decomposition
            method: Energy calculation method for "free" job type:
                - "pb/gbsa": Both PB and GB calculations
                - "gbsa": GB calculation only
            step_size: Frame sampling interval. Defaults to 10.

        Returns:
            File: Created input file for MMPBSA.py calculations.

        Raises:
            ValueError: If method is not specified for "free" job type or
                if invalid job_type/method is provided.

        Note:
            Input files saved in energy_inputfile/ subdirectory with
            appropriate parameter sets for different calculation types.
        """

        save_dir = os.path.join(self.save_dir, "energy_inputfile")
        os.makedirs(save_dir, exist_ok=True)

        if job_type == "entropy":
            self.inputfile = core.create_energy_inputfile(
                "nmode",
                startframe=start_frame,
                endframe=end_frame,
                step_size=step_size,
                save_dir=save_dir,
            )
        elif job_type == "free":
            if method is None:
                raise ValueError("Please specify method.")
            elif method == "pb/gbsa":
                self.inputfile = core.create_energy_inputfile(
                    "pb/gb",
                    startframe=start_frame,
                    endframe=end_frame,
                    step_size=step_size,
                    save_dir=save_dir,
                )
            elif method == "gbsa":
                self.inputfile = core.create_energy_inputfile(
                    "gb",
                    startframe=start_frame,
                    endframe=end_frame,
                    step_size=step_size,
                    save_dir=save_dir,
                )
            else:
                raise ValueError(f"Invalid method: {method}")
        elif job_type == "decomp":
            self.inputfile = core.create_energy_inputfile(
                "gb",
                startframe=start_frame,
                endframe=end_frame,
                step_size=step_size,
                decomp=True,
                save_dir=save_dir,
            )
        else:
            raise ValueError(f"Invalid job type: {job_type}")
        return self.inputfile

    def run_energy_calc(
        self,
        energy_decompose: str = None,
        cpu_num: int = None,
        input_file: File = None,
    ) -> None:
        """Execute MM/PBSA or MM/GBSA energy calculations.

        Runs AMBER's MMPBSA.py tool to calculate binding free energies,
        entropy contributions, or per-residue energy decomposition.

        Args:
            energy_decompose: Optional decomposition specification.
            cpu_num: Number of CPU cores to use. Defaults to system CPU_NUM.
            input_file: Custom input file to use. Defaults to self.inputfile.

        Raises:
            ValueError: If required input file or topology files are missing.

        Note:
            Requires all topology files (complex, receptor, ligand) and
            trajectory data to be loaded. Results saved in energy_results/
            subdirectory with detailed energy breakdown and statistics.
        """
        save_dir = os.path.join(self.save_dir, "energy_results")
        os.makedirs(save_dir, exist_ok=True)
        cpu_num = cpu_num if cpu_num is not None else CPU_NUM
        input_file = input_file if input_file is not None else self.inputfile

        if self.inputfile is None:
            raise ValueError("Please create inputfile first.")
        if self.top_file is None or self.traj_file is None:
            raise ValueError("Please load solvated complex top file first.")
        if not all([self.com_topfile, self.recep_topfile, self.ligand_topfile]):
            raise ValueError("Please load complex/receptor/ligand top file first.")

        core.run_energy_calculation(
            input_file=input_file,
            comsolvate_topfile=self.top_file,
            com_topfile=self.com_topfile,
            receptor_topfile=self.recep_topfile,
            ligand_topfile=self.ligand_topfile,
            traj_file=self.traj_file,
            save_dir=save_dir,
            energy_decompose=energy_decompose,
            cpu_num=cpu_num,
        )
        logger.info(f"Calculating normally finished.")
