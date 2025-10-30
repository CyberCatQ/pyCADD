"""High-level Gaussian quantum chemistry calculation interface.

This module provides the Gauss class which serves as a high-level interface for
running Gaussian quantum chemistry calculations including geometry optimization,
single-point energy calculations, and RESP/RESP2 charge calculations.
"""

import logging
import os
import tempfile

logger = logging.getLogger(__name__)

from pyCADD.Density.core import (
    _get_atom_lines,
    _get_mol2_lines_with_charge,
    generate_energy_input,
    generate_opt_input,
    merge_mol2_charge,
)
from pyCADD.utils.common import ChDir, File
from pyCADD.utils.tool import _find_execu, read_file, shell_run, write_file


class Gauss:
    """High-level interface for Gaussian quantum chemistry calculations.

    This class provides methods for running various types of Gaussian calculations
    including geometry optimization, single-point energy calculations, and RESP
    charge calculations. It automatically manages input file generation, job
    execution, and output file handling.

    Attributes:
        gauss (str): Path to the Gaussian executable (g16).

    Example:
        >>> gauss = Gauss()
        >>> opt_file = gauss.calc_opt("molecule.sdf", charge=0, multiplicity=1)
        >>> resp_file = gauss.calc_resp("molecule.sdf", charge=0)
    """

    def __init__(self) -> None:
        """Initializes the Gauss calculator by locating the Gaussian executable.

        Raises:
            FileNotFoundError: If the g16 executable is not found in the system PATH.
        """
        self.gauss = _find_execu("g16")
        if not self.gauss:
            raise FileNotFoundError("Gaussian executable 'g16' not found in system PATH.")

    @staticmethod
    def change_atom_type(mol2_block: str, atom_type: str = "gaff2", dumpto: str = None) -> str:
        with tempfile.TemporaryDirectory() as tmpdir:
            dumpto = dumpto or tmpdir
            with ChDir(dumpto):
                input_mol2_file = write_file("input.mol2", mol2_block)
                output_mol2_file = File("output.mol2", exist=False)
                shell_run(
                    f"antechamber -fi mol2 -i {input_mol2_file} -fo mol2 -o {output_mol2_file.file_path} -at {atom_type}"
                )
                return output_mol2_file.read()

    def run(self, gau_input: str, job_name: str = None, dumpto: str = None) -> str:
        """Executes a Gaussian calculation with the provided input.

        Runs a Gaussian calculation by writing the input to a file and executing
        the Gaussian program. The calculation can be run in a temporary directory
        or a specified output directory.

        Args:
            gau_input (str): Complete Gaussian input file content as a string.
            job_name (str, optional): Base name for input and output files. Defaults to "gauss_job".
            dumpto (str, optional): Directory path where calculation files should be saved.
                If None, uses a temporary directory that is automatically cleaned up.

        Returns:
            str: Content of the Gaussian output file.

        Raises:
            RuntimeError: If the Gaussian calculation fails.
            FileNotFoundError: If the Gaussian executable is not found.

        Note:
            When dumpto is None, all calculation files are created in a temporary
            directory and only the output content is returned. When dumpto is
            specified, all files are preserved in that directory.
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            work_dir = dumpto if dumpto else tmpdir
            os.makedirs(work_dir, exist_ok=True)
            with ChDir(work_dir):
                write_file(f"{job_name}.gjf", gau_input)
                job_name = job_name or "gauss_job"
                output_file = f"{job_name}.out"
                shell_run(f"{self.gauss} < {job_name}.gjf > {output_file}")
                return read_file(output_file)

    @staticmethod
    def convert_to_mol2_block(structure_file: str | File) -> str:
        """Converts a molecular structure file to MOL2 format.

        Uses OpenBabel to convert various molecular file formats (SDF, PDB, XYZ, etc.)
        to MOL2 format and returns the content as a string.

        Args:
            structure_file (str | File): Path to the input molecular structure file or File object.

        Returns:
            str: Complete MOL2 format content of the converted structure.

        Raises:
            RuntimeError: If OpenBabel conversion fails.
            FileNotFoundError: If the input structure file doesn't exist.

        Note:
            The conversion is performed in a temporary directory and only the
            MOL2 content is returned. Supported input formats include most
            common molecular file types recognized by OpenBabel.
        """
        structure_file = File(structure_file)
        with tempfile.TemporaryDirectory() as tmpdir:
            with ChDir(tmpdir):
                mol2_file = File(f"{structure_file.file_prefix}.mol2", exist=False)
                shell_run(
                    f"obabel -i{structure_file.file_ext} {structure_file.file_path} -o mol2 -O {mol2_file.file_path}"
                )
                return mol2_file.read()

    @staticmethod
    def get_resp_mol2_block(gout_file: str | File) -> str:
        """Extracts RESP charges from Gaussian output and generates MOL2 format.

        Uses AmberTools antechamber to extract RESP (Restrained Electrostatic
        Potential) charges from a Gaussian output file and create a MOL2 format
        structure with the calculated charges.

        Args:
            gout_file (str | File): Path to Gaussian output file (.out) containing ESP calculation
                results or File object.

        Returns:
            str: Complete MOL2 format content with RESP charges assigned to atoms.

        Raises:
            RuntimeError: If antechamber processing fails or if the Gaussian output
                doesn't contain the required ESP information.
            FileNotFoundError: If the Gaussian output file doesn't exist.

        Note:
            The input Gaussian calculation must have been performed with ESP
            calculation keywords (e.g., "pop=MK IOp(6/33=2,6/42=6)") to generate
            the electrostatic potential data required for RESP charge fitting.
        """
        gout_file = File(gout_file)
        with tempfile.TemporaryDirectory() as tmpdir:
            with ChDir(tmpdir):
                mol2_file = File(f"{gout_file.file_prefix}.mol2", exist=False)
                shell_run(
                    f"antechamber -i {gout_file.file_path} -fi gout"
                    f" -o {mol2_file.file_path} -fo mol2 -c resp"
                )
                return mol2_file.read()

    @staticmethod
    def get_resp2_mol2_block(
        gas_mol2_block: str, solvent_mol2_block: str, delta: float = 0.5
    ) -> str:
        """Generates RESP2 charges by interpolating between gas and solvent RESP charges.

        Calculates RESP2 charges using the formula:
        RESP2 = (1-δ) × RESP_gas + δ × RESP_solvent
        where δ (delta) is the interpolation parameter.

        Args:
            gas_mol2_block (str): MOL2 format string containing gas-phase RESP charges.
            solvent_mol2_block (str): MOL2 format string containing solvent-phase RESP charges.
            delta (float, optional): Interpolation parameter (0.0 = pure gas, 1.0 = pure solvent).
                Defaults to 0.5 for equal weighting.

        Returns:
            str: MOL2 format string with interpolated RESP2 charges.

        Raises:
            ValueError: If the number of atoms in gas and solvent MOL2 blocks don't match.
            IndexError: If either MOL2 block has incomplete atomic information.

        Note:
            Both input MOL2 blocks must represent the same molecular structure with
            identical atom ordering. The RESP2 method provides charges that balance
            gas-phase and solution-phase electrostatic properties.
        """
        gas_resp_charges = [float(line[8]) for line in _get_atom_lines(gas_mol2_block)]
        solvent_resp_charges = [float(line[8]) for line in _get_atom_lines(solvent_mol2_block)]
        if not len(gas_resp_charges) == len(solvent_resp_charges):
            raise ValueError("The number of atoms in gas and solvent mol2 blocks do not match.")
        adjusted_charges = [
            (1 - delta) * gas_charge + delta * solvent_charge
            for gas_charge, solvent_charge in zip(gas_resp_charges, solvent_resp_charges)
        ]
        adjusted_charges = [round(charge, 6) for charge in adjusted_charges]
        return _get_mol2_lines_with_charge(gas_mol2_block, adjusted_charges)

    def calc_opt(
        self,
        structure_file: str | File,
        charge: int,
        multiplicity: int = 1,
        dft: str = "B3LYP",
        basis_set: str = "6-31g*",
        solvent: str = "water",
        loose: bool = True,
        dispersion_correct: bool = False,
        td: bool = False,
        freq: bool = False,
        mem_use: str = "4GB",
        cpu_num: int = 4,
        save_dir: str = None,
    ) -> File:
        """Performs geometry optimization calculation using Gaussian.

        Executes a complete geometry optimization calculation with customizable
        DFT method, basis set, and optimization parameters. The optimized structure
        is saved to an output file.

        Args:
            structure_file (str | File): Path to the input molecular structure file or File object.
            charge (int): Molecular charge (required parameter).
            multiplicity (int, optional): Spin multiplicity. Defaults to 1.
            dft (str, optional): DFT functional name. Defaults to "B3LYP".
            basis_set (str, optional): Basis set specification. Defaults to "6-31g*".
            solvent (str, optional): Solvent for implicit solvation model. Defaults to "water".
            loose (bool, optional): Whether to use loose convergence criteria for faster optimization.
                Defaults to True.
            dispersion_correct (bool, optional): Whether to apply dispersion correction (GD3BJ).
                Defaults to False.
            td (bool, optional): Whether to perform time-dependent DFT optimization. Defaults to False.
            freq (bool, optional): Whether to calculate vibrational frequencies after optimization.
                Defaults to False.
            mem_use (str, optional): Memory allocation string. Defaults to "4GB".
            cpu_num (int, optional): Number of processors to use. Defaults to 4.
            save_dir (str, optional): Directory to save calculation result files. If None, saves in current directory.

        Returns:
            File: Generated Gaussian output file containing optimization results.

        Raises:
            RuntimeError: If the Gaussian optimization calculation fails.
            FileNotFoundError: If the input structure file doesn't exist.

        Note:
            The output file is named "{structure_prefix}_opt.out" and contains the
            complete optimization trajectory and final optimized geometry.
        """
        structure_file = File(structure_file)
        logger.info(f"Running geometry optimization for {structure_file.file_path}")
        logger.info(f"Charge = {charge}, Multiplicity = {multiplicity}, Solvent = {solvent}")
        opt_input = generate_opt_input(
            structure_file,
            charge=charge,
            multiplicity=multiplicity,
            dft=dft,
            basis_set=basis_set,
            solvent=solvent,
            loose=loose,
            dispersion_correct=dispersion_correct,
            td=td,
            freq=freq,
            mem_use=mem_use,
            cpu_num=cpu_num,
        )
        save_dir = save_dir or os.getcwd()
        os.makedirs(save_dir, exist_ok=True)
        opt_prefix = f"{structure_file.file_prefix}_opt"
        opt_output_file = os.path.join(save_dir, f"{opt_prefix}.out")
        opt_output = self.run(opt_input, "gauss_opt", dumpto=save_dir)
        logger.info(f"Optimization output saved: {opt_output_file}")
        return File(write_file(opt_output_file, opt_output))

    def calc_energy(
        self,
        structure_file: str | File,
        charge: int,
        multiplicity: int = 1,
        dft: str = "B3LYP",
        basis_set: str = "6-31g*",
        solvent: str = "water",
        dispersion_correct: bool = False,
        td: bool = False,
        esp_calculate: bool = False,
        mem_use: str = "4GB",
        cpu_num: int = 4,
        save_dir: str = None,
    ) -> File:
        """Performs single-point energy calculation using Gaussian.

        Calculates the electronic energy of a molecular structure without geometry
        optimization. Can include electrostatic potential calculation for RESP
        charge fitting and time-dependent DFT for excited states.

        Args:
            structure_file (str | File): Path to the input molecular structure file or File object.
            charge (int): Molecular charge (required parameter).
            multiplicity (int, optional): Spin multiplicity. Defaults to 1.
            dft (str, optional): DFT functional name. Defaults to "B3LYP".
            basis_set (str, optional): Basis set specification. Defaults to "6-31g*".
            solvent (str, optional): Solvent for implicit solvation model. Defaults to "water".
                Set to None for gas-phase calculation.
            dispersion_correct (bool, optional): Whether to apply dispersion correction (GD3BJ).
                Defaults to False.
            td (bool, optional): Whether to perform time-dependent DFT calculation. Defaults to False.
            esp_calculate (bool, optional): Whether to calculate electrostatic potential for RESP
                charge fitting. Defaults to False.
            mem_use (str, optional): Memory allocation string. Defaults to "4GB".
            cpu_num (int, optional): Number of processors to use. Defaults to 4.
            save_dir (str, optional): Directory to save calculation result files. If None, saves in current directory.

        Returns:
            File: Path to the generated Gaussian output file containing energy results.

        Raises:
            RuntimeError: If the Gaussian energy calculation fails.
            FileNotFoundError: If the input structure file doesn't exist.

        Note:
            Output file naming follows the pattern "{structure_prefix}_{solvent}_energy.out"
            or "{structure_prefix}_gas_energy.out" for gas-phase calculations.
            When esp_calculate=True, the output contains ESP data suitable for RESP fitting.
        """
        structure_file = File(structure_file)
        logger.info(f"Running single point energy calculation for {structure_file.file_path}")
        logger.info(
            f"Charge = {charge}, Multiplicity = {multiplicity}, Solvent = {solvent if solvent else 'None(Gas)'}"
        )
        energy_input = generate_energy_input(
            structure_file,
            charge=charge,
            multiplicity=multiplicity,
            dft=dft,
            basis_set=basis_set,
            solvent=solvent,
            dispersion_correct=dispersion_correct,
            td=td,
            esp_calculate=esp_calculate,
            mem_use=mem_use,
            cpu_num=cpu_num,
        )
        save_dir = save_dir or os.getcwd()
        os.makedirs(save_dir, exist_ok=True)
        ene_prefix = f"{structure_file.file_prefix}"
        ene_prefix += f"_{solvent}" if solvent else "_gas"
        ene_prefix += "_energy"
        ene_result_file_path = os.path.join(save_dir, f"{ene_prefix}.out")
        energy_output = self.run(energy_input, "gauss_energy", dumpto=save_dir)
        logger.info(f"Energy calculation output saved: {ene_result_file_path}")
        return File(write_file(ene_result_file_path, energy_output))

    def calc_resp(
        self,
        structure_file: str | File,
        charge: int,
        multiplicity: int = 1,
        dft: str = "B3LYP",
        basis_set: str = "6-31g*",
        solvent: str = "water",
        mem_use: str = "4GB",
        cpu_num: int = 4,
        atom_type: str = "gaff2",
        save_dir: str = None,
    ) -> File:
        """Performs complete RESP charge calculation workflow.

        Executes a full RESP (Restrained Electrostatic Potential) charge calculation
        workflow including geometry optimization and single-point energy calculation
        with ESP analysis. The resulting charges are merged back into the original
        molecular structure format.

        Args:
            structure_file (str | File): Path to the input molecular structure file or File object.
            charge (int): Molecular charge (required parameter).
            multiplicity (int, optional): Spin multiplicity. Defaults to 1.
            dft (str, optional): DFT functional name. Defaults to "B3LYP".
            basis_set (str, optional): Basis set specification. Defaults to "6-31g*".
            solvent (str, optional): Solvent for implicit solvation model. Defaults to "water".
            mem_use (str, optional): Memory allocation string. Defaults to "4GB".
            cpu_num (int, optional): Number of processors to use. Defaults to 4.
            atom_type (str, optional): Atom type to assign using antechamber. Defaults to "gaff2". None for no change.
            save_dir (str, optional): Directory to save calculation result files. If None, saves in current directory.

        Returns:
            File: Path to the generated MOL2 file with RESP charges assigned to atoms.

        Raises:
            RuntimeError: If any step of the RESP calculation workflow fails.
            FileNotFoundError: If the input structure file doesn't exist.

        Note:
            The complete workflow includes:
            1. Geometry optimization with loose convergence
            2. Single-point energy calculation with ESP analysis
            3. RESP charge fitting using antechamber
            4. Merging charges back into original molecular structure
            Output file is named "{structure_prefix}_resp.mol2".
        """
        structure_file = File(structure_file)
        logger.info(f"Running RESP calculation for {structure_file.file_path}")
        save_dir = save_dir or os.getcwd()
        os.makedirs(save_dir, exist_ok=True)
        opt_output_file = self.calc_opt(
            structure_file,
            charge=charge,
            multiplicity=multiplicity,
            dft=dft,
            basis_set=basis_set,
            loose=True,
            dispersion_correct=False,
            mem_use=mem_use,
            cpu_num=cpu_num,
            save_dir=save_dir,
        )
        gas_resp_output_file = self.calc_energy(
            opt_output_file,
            charge=charge,
            multiplicity=multiplicity,
            dft=dft,
            basis_set=basis_set,
            solvent=solvent,
            dispersion_correct=False,
            td=False,
            mem_use=mem_use,
            cpu_num=cpu_num,
            esp_calculate=True,
            save_dir=save_dir,
        )
        ori_mol2 = self.convert_to_mol2_block(structure_file)
        resp_mol2 = self.get_resp_mol2_block(gas_resp_output_file)
        merged_mol2 = merge_mol2_charge(ori_mol2, resp_mol2)
        merged_mol2_file_path = os.path.join(save_dir, f"{structure_file.file_prefix}_resp.mol2")
        if atom_type is not None:
            merged_mol2 = self.change_atom_type(merged_mol2, atom_type=atom_type)
        logger.info(f"RESP output saved: {merged_mol2_file_path}")
        return File(write_file(merged_mol2_file_path, merged_mol2))

    def calc_resp2(
        self,
        structure_file: str | File,
        charge: int,
        multiplicity: int = 1,
        dft: str = "B3LYP",
        basis_set: str = "6-31g*",
        solvent: str = "water",
        mem_use: str = "4GB",
        cpu_num: int = 4,
        delta: float = 0.5,
        atom_type: str = "gaff2",
        save_dir: str = None,
    ) -> File:
        """Performs complete RESP2 charge calculation workflow.

        Executes a full RESP2 charge calculation workflow that combines gas-phase
        and solvent-phase RESP charges using interpolation. This method provides
        charges that balance gas-phase and solution-phase electrostatic properties.

        Args:
            structure_file (str | File): Path to the input molecular structure file or File object.
            charge (int): Molecular charge (required parameter).
            multiplicity (int, optional): Spin multiplicity. Defaults to 1.
            dft (str, optional): DFT functional name. Defaults to "B3LYP".
            basis_set (str, optional): Basis set specification. Defaults to "6-31g*".
            solvent (str, optional): Solvent for implicit solvation model. Defaults to "water".
            mem_use (str, optional): Memory allocation string. Defaults to "4GB".
            cpu_num (int, optional): Number of processors to use. Defaults to 4.
            delta (float, optional): Interpolation parameter for RESP2 charge calculation.
                0.0 = pure gas-phase, 1.0 = pure solvent-phase. Defaults to 0.5.
            atom_type (str, optional): Atom type to assign using antechamber. Defaults to "gaff2". None for no change.
            save_dir (str, optional): Directory to save calculation result files. If None, saves in current directory.

        Returns:
            File: Path to the generated MOL2 file with RESP2 charges assigned to atoms.

        Raises:
            RuntimeError: If any step of the RESP2 calculation workflow fails.
            FileNotFoundError: If the input structure file doesn't exist.
            ValueError: If gas and solvent calculations produce inconsistent atom counts.

        Note:
            The complete workflow includes:
            1. Geometry optimization with loose convergence
            2. Gas-phase single-point energy calculation with ESP analysis
            3. Solvent-phase single-point energy calculation with ESP analysis
            4. RESP charge fitting for both phases using antechamber
            5. Interpolation of charges using RESP2 formula: (1-δ)*RESP_gas + δ*RESP_solvent
            Output files include individual gas/solvent RESP files and final RESP2 file.
        """
        structure_file = File(structure_file)
        logger.info(f"Running RESP2({delta}) calculation for {structure_file.file_path}")
        save_dir = save_dir or os.getcwd()
        os.makedirs(save_dir, exist_ok=True)
        original_mol2_block = self.convert_to_mol2_block(structure_file)
        opt_output_file = self.calc_opt(
            structure_file=structure_file,
            charge=charge,
            multiplicity=multiplicity,
            dft=dft,
            basis_set=basis_set,
            loose=True,
            dispersion_correct=False,
            mem_use=mem_use,
            cpu_num=cpu_num,
            save_dir=save_dir,
        )
        gas_resp_output_file = self.calc_energy(
            structure_file=opt_output_file,
            charge=charge,
            multiplicity=multiplicity,
            dft=dft,
            basis_set=basis_set,
            solvent=None,  # No solvent for gas phase
            dispersion_correct=False,
            td=False,
            mem_use=mem_use,
            cpu_num=cpu_num,
            esp_calculate=True,
            save_dir=save_dir,
        )
        gas_mol2_block = merge_mol2_charge(
            self.get_resp_mol2_block(gas_resp_output_file), original_mol2_block
        )
        gas_resp_mol2_file_path = os.path.join(
            save_dir, f"{structure_file.file_prefix}_gas_resp.mol2"
        )
        write_file(gas_resp_mol2_file_path, gas_mol2_block)
        logger.info(f"Gas-phase energy calculation output saved: {gas_resp_mol2_file_path}")
        solvent_resp_output_file = self.calc_energy(
            structure_file=opt_output_file,
            charge=charge,
            multiplicity=multiplicity,
            dft=dft,
            basis_set=basis_set,
            solvent=solvent,
            dispersion_correct=False,
            td=False,
            mem_use=mem_use,
            cpu_num=cpu_num,
            esp_calculate=True,
            save_dir=save_dir,
        )
        solvent_mol2_block = merge_mol2_charge(
            self.get_resp_mol2_block(solvent_resp_output_file), original_mol2_block
        )
        solvent_resp_mol2_file_path = os.path.join(
            save_dir, f"{structure_file.file_prefix}_solvent_resp.mol2"
        )
        write_file(solvent_resp_mol2_file_path, solvent_mol2_block)
        resp2_mol2_file_path = os.path.join(save_dir, f"{structure_file.file_prefix}_resp2.mol2")
        logger.info(f"Solvent-phase energy calculation output saved: {solvent_resp_mol2_file_path}")
        resp2_mol2_block = self.get_resp2_mol2_block(
            gas_mol2_block=gas_mol2_block,
            solvent_mol2_block=solvent_mol2_block,
            delta=delta,
        )
        if atom_type is not None:
            resp2_mol2_block = self.change_atom_type(resp2_mol2_block, atom_type=atom_type)
        result_file = File(
            write_file(
                resp2_mol2_file_path,
                resp2_mol2_block,
            )
        )
        logger.info(f"RESP2 output saved: {resp2_mol2_file_path}")
        return result_file
