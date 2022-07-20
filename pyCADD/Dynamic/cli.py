import os
import warnings

import click

warnings.filterwarnings("ignore")

@click.command()
@click.argument('protein_file', type=click.Path(exists=True))
@click.argument('molecule_file', type=click.Path(exists=True))
@click.option('--charge', '-c', default=0, show_default=True, type=int, help='Charge of molecule.')
@click.option('--multiplicity', '-m', default=1, show_default=True, type=int, help='Multiplicity of molecule.')
@click.option('--solvent', '-s', default='water', show_default=True, type=str, help='Solvent used in molecule preparation.')
@click.option('--prefix', '-p', default=None, type=str, help='Prefix of output files.')
@click.option('--parallel', '-n', default=os.cpu_count(), show_default=True, type=int, help='Number of parallel processes.')
@click.option('--with-gpu', '-g', default=0, show_default=True, type=int, help='Specify GPU device code used in simulation. Default to 0')
def main(protein_file, molecule_file, charge, multiplicity, solvent, prefix, parallel, with_gpu):
    '''
    Prepare necessary files and run molecular dynamics simulation.\n
    protein_file : Specify protein file (PDB format) path for MD.\n
    molecule_file : Specify molecule file (PDB format) path for MD.
    '''
    from pyCADD.Dynamic import Processor, Simulator
    processor = Processor()
    processor.protein_prepare(protein_file)
    processor.molecule_prepare(
        molecule_file, charge, multiplicity, parallel, solvent)
    prefix = os.path.basename(os.getcwd()) if prefix is None else prefix
    processor.leap_prepare(prefix)
    simulator = Simulator(processor)
    simulator.creat_input_file()
    simulator.run_simulation(with_gpu)


@click.command()
@click.option('-y', type=click.Path(exists=True), help='Trajectory file path.', prompt='Please specify trajectory file path.')
@click.option('-sp', type=click.Path(exists=True), help='Solvated complex topology file path.', prompt='Please specify solvated complex topology file path:')
@click.option('-lp', type=click.Path(exists=True), help='Ligand topology file path.', prompt='Please specify ligand topology file path:')
@click.option('-rp', type=click.Path(exists=True), help='Receptor topology file path.', prompt='Please specify receptor topology file path:')
@click.option('-cp', type=click.Path(exists=True), help='Complex topology file path.', prompt='Please specify complex topology file path:')
@click.option('-ro', type=click.Path(exists=True), help='Molecular dynamics output file path.', prompt='Please specify molecular dynamics output file path:')
@click.option('--no-hbond', '-nh', is_flag=True, help='Disable calculating and tracing of hydrogen bonds.')
@click.option('--no-rmsd', '-nd', is_flag=True, help='Disable calculating of RMSD.')
@click.option('--no-rmsf', '-nf', is_flag=True, help='Disable calculating of RMSF.')
@click.option('--no-extract', '-ne', is_flag=True, help='Disable extracting of lowest energy structure.')
@click.option(
    '--decomP', '-d', 
    type=(int, int, int), 
    help='Performing MM-GBSA energy decomposition from START_FRAME to END_FRAME with STEP_SIZE.', 
    prompt='Please specify decomposition START_FRAME, END_FRAME, STEP_SIZE:'
    )
@click.option(
    '--nmode', '-e', 
    type=(int, int, int), 
    help='Performing entropy calculation with normal mode(nmode) from START_FRAME to END_FRAME with STEP_SIZE.', 
    prompt='Please specify entropy calculation START_FRAME, END_FRAME, STEP_SIZE:'
    )
@click.option('--parallel', '-n', default=os.cpu_count(), show_default=True, type=int, help='Number of parallel processes used in energy calculation.')
def analysis(y, sp, lp, rp, cp, ro, no_hbond, no_rmsd, no_rmsf, no_extract, decomp=None, nmode=None, parallel=None):
    '''
    Post-process for molecular dynamics simulation.\n
    Workflow:\n
        1. Calculate RMSD\n
        2. Calculate RMSF\n
        3. Extract lowest energy structure\n
        4. Calculate and trace hydrogen bonds\n
        [Optional]\n
        * Perform MM-GBSA energy decomposition\n
        * Calculate entropy with normal mode\n
    '''
    from pyCADD.Dynamic import Analyzer
    parallel = os.cpu_count() if parallel is None else parallel
    analyzer = Analyzer(
        traj_file_path=y,
        comsolvated_topfile_path=sp,
        com_topfile_path=cp,
        ligand_topfile_path=lp,
        receptor_topfile_path=rp,
        mdout_file_path=ro
        )

    if not no_rmsd:
        analyzer.calc_rmsd()
    if not no_rmsf:
        analyzer.calc_rmsf()
    if not no_hbond:
        analyzer.calc_hbond()
    if not no_extract:
        analyzer.extract_lowest_energy_st()
    if decomp is not None:
        start_frame, end_frame, step_size = decomp
        analyzer.creat_energy_inputfile(start_frame=start_frame, end_frame=end_frame, interval=step_size, job_type='decomp')
        analyzer.run_energy_calc(cpu_num=parallel)
    if nmode is not None:
        start_frame, end_frame, step_size = nmode
        analyzer.creat_energy_inputfile(start_frame=start_frame, end_frame=end_frame, interval=step_size, job_type='entropy')
        analyzer.run_energy_calc(cpu_num=parallel)
