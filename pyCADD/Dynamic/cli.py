import os
import warnings

import click

warnings.filterwarnings("ignore")

from pyCADD.Dynamic import Processor, Simulator, Analyzer


def md_prepare(
    protein_file: str,
    molecule_file: str,
    charge: int = 0,
    multiplicity: int = 1,
    solvent: str = "water",
    parallel: int = 4,
    time: float = 100.0,
    dft: str = "B3LYP",
    basis_set: str = "6-31g*",
    charge_method: str = "resp",
    keep_water: bool = False,
    box_size: float = 10.0,
    overwrite: bool = False,
):
    apo = molecule_file is None
    processor = Processor(apo=apo)
    processor.molecule_prepare(
        molecule_file_path=molecule_file,
        dft=dft,
        basis_set=basis_set,
        charge=charge,
        multiplicity=multiplicity,
        cpu_num=parallel,
        charge_method=charge_method,
        solvent=solvent,
        overwrite=overwrite,
    )
    processor.protein_prepare(
        protein_file_path=protein_file, keep_water=keep_water, overwrite=overwrite
    )
    processor.leap_prepare(box_size=box_size)
    step_num = int(time * 1000 / 0.002)
    water_resnum = processor.get_water_resnum()
    water_resnum_start = int(water_resnum[0])
    water_resnum_end = int(water_resnum[-1])

    processor.add_minimize_process(
        process_name="minimize_complex",
        restraint=True,
        restraint_mask=f"':1-{water_resnum_start-1}'",
        use_gpu=False,
        cpu_num=parallel,
    )
    processor.add_minimize_process(
        process_name="minimize_solvent",
        restraint=True,
        restraint_mask=f"':{water_resnum_start}-{water_resnum_end}'",
        use_gpu=False,
        cpu_num=parallel,
    )
    processor.add_minimize_process(process_name="minimize_all", use_gpu=False, cpu_num=parallel)
    processor.add_heat_process()
    restraintmask = "'!(:WAT,Na+,Cl-,K+,K) & !@H= & !@H'"
    for rest_wt in [4.0, 3.5, 3.0, 2.5, 2.0, 1.0, 0]:
        processor.add_npt_process(
            total_step=5000,
            process_name=f"eq_npt_reswt{rest_wt}",
            restraint_wt=rest_wt,
            restraintmask=restraintmask,
        )
    processor.add_npt_process(total_step=500000, process_name="eq_npt")
    processor.add_npt_process(total_step=step_num, is_production=True, process_name="production")
    return processor


def md_simulate(top_file: str, inpcrd_file: str, with_gpu: int = 0):
    processor = Processor()
    processor.set_prepared_file(top_file, "com_top")
    processor.set_prepared_file(inpcrd_file, "com_crd")
    processor.add_process("input_file/minimize_complex.in", "minimize_complex", "minimize")
    processor.add_process("input_file/minimize_solvent.in", "minimize_solvent", "minimize")
    processor.add_process("input_file/minimize_all.in", "minimize_all", "minimize")
    processor.add_process("input_file/heat.in", "heat")
    for rest_wt in [4.0, 3.5, 3.0, 2.5, 2.0, 1.0, 0]:
        processor.add_process(
            f"input_file/eq_npt_reswt{rest_wt}.in", f"eq_npt_reswt{rest_wt}", "npt"
        )
    processor.add_process("input_file/eq_npt.in", "eq_npt", "npt")
    processor.add_process("input_file/production.in", "production", "npt")
    simulator = Simulator(processor)
    simulator.run_simulation(with_gpu)
    return simulator


def md_workflow(
    protein_file: str,
    molecule_file: str,
    charge: int = 0,
    multiplicity: int = 1,
    solvent: str = "water",
    parallel: int = 1,
    time: float = 1.0,
    dft: str = "B3LYP",
    basis_set: str = "6-31G*",
    charge_method: str = "RESP",
    keep_water: bool = False,
    box_size: float = 10.0,
    overwrite: bool = False,
    with_gpu: int = 0,
    analysis: bool = True,
):
    processor = md_prepare(
        protein_file,
        molecule_file,
        charge,
        multiplicity,
        solvent,
        parallel,
        time,
        dft,
        basis_set,
        charge_method,
        keep_water,
        box_size,
        overwrite,
    )
    simulator = Simulator(processor)
    simulator.run_simulation(cuda_device=with_gpu)
    if analysis:
        params = {
            "traj_file": simulator.traj_file,
            "comsolvated_topfile": simulator.comsolvate_topfile,
            "com_topfile": simulator.com_topfile,
            "receptor_topfile": simulator.pro_topfile,
            "ligand_topfile": simulator.lig_topfile,
            "mdout_file": simulator.mdout_file,
        }
        analyzer = Analyzer(save_dir="md_analysis", **params)
        analyzer.calc_rmsd()
        analyzer.calc_rmsf()
        analyzer.calc_hbond()
        analyzer.create_energy_inputfile(
            start_frame=0, end_frame=analyzer.traj.shape[0], step_size=10, job_type="decomp"
        )
        analyzer.run_energy_calc(cpu_num=parallel)


@click.group()
def main():
    """
    Command line interface for pyCADD-Dynamic.
    """
    pass


@main.command(short_help="Automatically preparing required files for MD only.")
@click.argument("protein_file", type=click.Path(exists=True))
@click.argument("molecule_file", type=click.Path(exists=True), required=False)
@click.option(
    "--charge", "-c", default=0, show_default=True, type=int, help="Net charge of molecule."
)
@click.option(
    "--multiplicity", "-m", default=1, show_default=True, type=int, help="Multiplicity of molecule."
)
@click.option(
    "--solvent",
    "-s",
    default="water",
    show_default=True,
    type=str,
    help="Solvent used in molecule preparation (RESP/RESP2).",
)
@click.option(
    "--parallel", "-n", default=4, show_default=True, type=int, help="Number of parallel processes."
)
@click.option(
    "--with-gpu",
    "-g",
    default=0,
    show_default=True,
    type=int,
    help="Specify GPU device code used in simulation. Default to 0",
)
@click.option(
    "--time",
    "-t",
    default=100,
    show_default=True,
    type=int,
    help="Total time(ns) of simulation. Default to 100 ns.",
)
@click.option(
    "--dft", "-d", default="B3LYP", show_default=True, type=str, help="DFT functional to use."
)
@click.option(
    "--basis-set",
    "-bs",
    default="6-31G*",
    show_default=True,
    type=str,
    help="Basis set to use in DFT calculation.",
)
@click.option(
    "--charge-method",
    "-cm",
    default="resp",
    show_default=True,
    type=click.Choice(["resp", "resp2", "bcc"], case_sensitive=False),
    help="Charge calculation method for molecule.",
)
@click.option(
    "--keep-water", "-w", is_flag=True, help="Keep water molecules in the original protein file."
)
@click.option(
    "--box-size",
    "-b",
    default=10.0,
    show_default=True,
    type=float,
    help="TIP3P Water Box size of simulation. Default to 10 Angstrom.",
)
@click.option("--overwrite", "-O", is_flag=True, help="Overwrite existing files.")
def prepare(
    protein_file,
    molecule_file,
    charge,
    multiplicity,
    solvent,
    parallel,
    time,
    dft,
    basis_set,
    charge_method,
    keep_water,
    box_size,
    overwrite,
):
    """
    Prepare required files for MD.\n
    protein_file : Specify protein file (PDB format) path for MD.\n
    molecule_file : Specify molecule file (PDB format) path for MD. If not specified, apo system(protein only) will be simulated.\n
    """
    md_prepare(
        protein_file=protein_file,
        molecule_file=molecule_file,
        charge=charge,
        multiplicity=multiplicity,
        solvent=solvent,
        parallel=parallel,
        time=time,
        dft=dft,
        basis_set=basis_set,
        charge_method=charge_method,
        keep_water=keep_water,
        box_size=box_size,
        overwrite=overwrite,
    )


@main.command(short_help="Running molecular dynamics simulation with prepared files.")
@click.argument("top_file", type=click.Path(exists=True))
@click.argument("inpcrd_file", type=click.Path(exists=True))
@click.option(
    "--with-gpu",
    "-g",
    default=0,
    show_default=True,
    type=int,
    help="Specify GPU device code used in simulation. Default to 0",
)
def simulate(top_file, inpcrd_file, with_gpu):
    """
    Run molecular dynamics simulation.\nOnly supports running workflows automatically generated by prepare.
    top_file : Specify solvated complex topology file (TOP format) path for MD.\n
    inpcrd_file : Specify solvated complex input coordinate file (INPCRD format) path for MD.
    """
    md_simulate(top_file=top_file, inpcrd_file=inpcrd_file, with_gpu=with_gpu)


@main.command(short_help="Automatic MD workflow (including preparation and simulation).")
@click.argument("protein_file", type=click.Path(exists=True))
@click.argument("molecule_file", type=click.Path(exists=True), required=False)
@click.option(
    "--charge", "-c", default=0, show_default=True, type=int, help="Net charge of molecule."
)
@click.option(
    "--multiplicity", "-m", default=1, show_default=True, type=int, help="Multiplicity of molecule."
)
@click.option(
    "--solvent",
    "-s",
    default="water",
    show_default=True,
    type=str,
    help="Solvent used in molecule preparation (RESP/RESP2).",
)
@click.option(
    "--parallel", "-n", default=4, show_default=True, type=int, help="Number of parallel processes."
)
@click.option(
    "--with-gpu",
    "-g",
    default=0,
    show_default=True,
    type=int,
    help="Specify GPU device code used in simulation. Default to 0",
)
@click.option(
    "--time",
    "-t",
    default=100,
    show_default=True,
    type=int,
    help="Total time(ns) of simulation. Default to 100 ns.",
)
@click.option(
    "--dft", "-d", default="B3LYP", show_default=True, type=str, help="DFT functional to use."
)
@click.option(
    "--basis-set",
    "-bs",
    default="6-31G*",
    show_default=True,
    type=str,
    help="Basis set to use in DFT calculation.",
)
@click.option(
    "--charge-method",
    "-cm",
    default="resp",
    show_default=True,
    type=click.Choice(["resp", "resp2", "bcc"], case_sensitive=False),
    help="Charge calculation method for molecule.",
)
@click.option(
    "--keep-water", "-w", is_flag=True, help="Keep water molecules in the original protein file."
)
@click.option(
    "--box-size",
    "-b",
    default=10,
    show_default=True,
    type=float,
    help="TIP3P Water Box size of simulation. Default to 10 Angstrom.",
)
@click.option("--analysis", "-a", is_flag=True, help="Perform analysis after simulation done.")
@click.option("--overwrite", "-O", is_flag=True, help="Overwrite existing files.")
def auto(
    protein_file,
    molecule_file,
    charge,
    multiplicity,
    solvent,
    parallel,
    with_gpu,
    time,
    dft,
    basis_set,
    charge_method,
    keep_water,
    box_size,
    analysis,
    overwrite,
):
    """
    Prepare required files and run molecular dynamics simulation.\n
    protein_file : Specify protein file (PDB format) path for MD.\n
    molecule_file : Specify molecule file (PDB format) path for MD. If not specified, apo system(protein only) will be simulated.\n
    """
    md_workflow(
        protein_file=protein_file,
        molecule_file=molecule_file,
        charge=charge,
        multiplicity=multiplicity,
        solvent=solvent,
        parallel=parallel,
        time=time,
        dft=dft,
        basis_set=basis_set,
        charge_method=charge_method,
        keep_water=keep_water,
        box_size=box_size,
        overwrite=overwrite,
        with_gpu=with_gpu,
        analysis=analysis,
    )


@main.command(short_help="Simple post-analysis for finished MD simulation.")
@click.option(
    "-y",
    type=click.Path(exists=True),
    help="Trajectory file path.",
    prompt="Please specify trajectory file path",
)
@click.option(
    "-sp",
    type=click.Path(exists=True),
    help="Solvated complex topology file path.",
    prompt="Please specify solvated complex topology file path",
)
@click.option(
    "-lp",
    type=click.Path(exists=True),
    help="Ligand topology file path. If not specified, apo system analysis will be performed.",
)
@click.option(
    "-rp",
    type=click.Path(exists=True),
    help="Receptor topology file path.",
    prompt="Please specify receptor topology file path",
)
@click.option(
    "-cp",
    type=click.Path(exists=True),
    help="Complex topology file path.",
    prompt="Please specify complex topology file path",
)
@click.option(
    "-ro",
    type=click.Path(exists=True),
    help="Molecular dynamics output file path.",
    prompt="Please specify molecular dynamics output file path",
)
@click.option(
    "--no-hbond", "-nh", is_flag=True, help="Disable calculating and tracing of hydrogen bonds."
)
@click.option("--no-rmsd", "-nd", is_flag=True, help="Disable calculating of RMSD.")
@click.option("--no-rmsf", "-nf", is_flag=True, help="Disable calculating of RMSF.")
@click.option(
    "--decomp",
    "-d",
    nargs=3,
    help="Performing MM-GBSA energy decomposition from START_FRAME(INT1) to END_FRAME(INT2) with STEP_SIZE(INT3).",
    default=None,
    type=int,
)
@click.option(
    "--nmode",
    "-e",
    nargs=3,
    help="Performing entropy calculation with normal mode(nmode) from START_FRAME(INT1) to END_FRAME(INT2) with STEP_SIZE(INT3).",
    default=None,
    type=int,
)
@click.option(
    "--parallel",
    "-n",
    default=4,
    show_default=True,
    type=int,
    help="Number of parallel processes used in energy calculation.",
)
def analysis(y, sp, lp, rp, cp, ro, no_hbond, no_rmsd, no_rmsf, decomp, nmode, parallel):
    """
    Post-process for molecular dynamics simulation.\n
    Workflow:\n
        1. Calculate RMSD\n
        2. Calculate RMSF\n
        3. Calculate and trace distance, angle and lifetime of hydrogen bonds\n
        [Optional]\n
        * Perform MM-GBSA energy decomposition\n
        * Calculate entropy with normal mode\n
    """
    analyzer = Analyzer(
        traj_file_path=y,
        comsolvated_topfile_path=sp,
        com_topfile_path=cp,
        ligand_topfile_path=lp,
        receptor_topfile_path=rp,
        mdout_file_path=ro,
        save_dir="md_analysis",
    )

    if decomp is not None:
        try:
            decomp_start_fm, decomp_end_fm, decomp_step_size = decomp
        except ValueError:
            raise ValueError(
                "Invalid input. Energy decomposition needs 3 integers for START_FRAME, END_FRAME and STEP_SIZE."
            )
        if decomp_start_fm > decomp_end_fm:
            raise ValueError("Invalid input. Energy decomposition can not be performed.")
        decomp = True

    if nmode is not None:
        try:
            nmode_start_fm, nmode_end_fm, nmode_step_size = nmode
        except ValueError:
            raise ValueError(
                "Invalid input. Entropy calculation needs 3 integers for START_FRAME, END_FRAME and STEP_SIZE."
            )
        if nmode_start_fm > nmode_end_fm:
            raise ValueError("Invalid input. Entropy calculation can not be performed.")
        nmode = True

    if not no_rmsd:
        analyzer.calc_rmsd()
    if not no_rmsf:
        analyzer.calc_rmsf()
    if not no_hbond:
        analyzer.calc_hbond()
    if decomp:
        analyzer.create_energy_inputfile(
            start_frame=decomp_start_fm,
            end_frame=decomp_end_fm,
            step_size=decomp_step_size,
            job_type="decomp",
        )
        analyzer.run_energy_calc(cpu_num=parallel)
    if nmode:
        analyzer.create_energy_inputfile(
            start_frame=nmode_start_fm,
            end_frame=nmode_end_fm,
            step_size=nmode_step_size,
            job_type="entropy",
        )
        analyzer.run_energy_calc(cpu_num=parallel)
