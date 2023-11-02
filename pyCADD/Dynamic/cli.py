import os
import warnings

import click

warnings.filterwarnings("ignore")

@click.group()
def main():
    '''
    Command line interface for pyCADD-Dynamic.
    '''
    pass

@main.command(short_help='Automatic MD workflow (including preparation and simulation).')
@click.argument('protein_file', type=click.Path(exists=True))
@click.argument('molecule_file', type=click.Path(exists=True), required=False)
@click.option('--charge', '-c', default=0, show_default=True, type=int, help='Charge of molecule.')
@click.option('--multiplicity', '-m', default=1, show_default=True, type=int, help='Multiplicity of molecule.')
@click.option('--solvent', '-s', default='water', show_default=True, type=str, help='Solvent used in molecule preparation.')
@click.option('--prefix', '-p', default=None, type=str, help='Prefix of LEaP output files.')
@click.option('--parallel', '-n', default=os.cpu_count(), show_default=True, type=int, help='Number of parallel processes.')
@click.option('--with-gpu', '-g', default=0, show_default=True, type=int, help='Specify GPU device code used in simulation. Default to 0')
@click.option('--time', '-t', default=100, show_default=True, type=int, help='Total time(ns) of simulation. Default to 100 ns.')
@click.option('-bcc', is_flag=True, help='Use existing BCC charges instead of RESP.')
@click.option('--keep-cood', '-k', is_flag=True, help='Keep atoms coordinates from the input molecular file after calculating resp charge, instead of using the gaussian optimized output.')
@click.option('--keep-water', '-w', is_flag=True, help='Keep water molecules in the original protein file.')
@click.option('--box-size', '-b', default=12, show_default=True, type=float, help='TIP3P Water Box size of simulation. Default to 12 Angstrom.')
@click.option('--overwrite', '-O', is_flag=True, help='Overwrite existing files.')
def auto(protein_file, molecule_file, charge, multiplicity, solvent, prefix, parallel, with_gpu, time, bcc, keep_cood, keep_water, box_size, overwrite):
    '''
    Prepare required files and run molecular dynamics simulation.\n
    protein_file : Specify protein file (PDB format) path for MD.\n
    molecule_file : Specify molecule file (PDB format) path for MD. If not specified, apo system(protein only) will be simulated.\n
    '''
    from pyCADD.Dynamic import Processor, Simulator
    # Prepare files
    processor = Processor() if molecule_file is not None else Processor(apo=True)
    if molecule_file is not None:
        if not bcc:
            processor.molecule_prepare(
                molecule_file, charge, multiplicity, parallel, solvent, overwrite=overwrite, method='resp', keep_origin_cood=keep_cood)
        else:
            processor.molecule_prepare(
                molecule_file, charge, multiplicity, parallel, solvent, overwrite=overwrite, method='bcc')
            
    processor.protein_prepare(protein_file, keep_water=keep_water)
    prefix = os.path.basename(os.getcwd()) if prefix is None else prefix
    processor.leap_prepare(prefix=prefix, box_size=box_size)
    step_num = int(time * 1000 / 0.002)
    
    water_resnum = processor.get_water_resnum()
    water_resnum_start = int(water_resnum[0])
    water_resnum_end = int(water_resnum[-1])
    
    processor.add_minimize_process(process_name='stepA', restraint=True, restraint_mask=f"':1-{water_resnum_start-1}'")
    processor.add_minimize_process(process_name='stepB', restraint=True, restraint_mask=f"':{water_resnum_start}-{water_resnum_end}'")
    processor.add_minimize_process(process_name='stepC')
    processor.add_heat_process()
    restraintmask = "'!(:WAT,Na+,Cl-,K+,K) & !@H= & !@H'"
    for rest_wt in [4.0, 3.5, 3.0, 2.5, 2.0, 1.0, 0]:
        processor.add_npt_process(total_step=5000, process_name=f'eq_npt_reswt{rest_wt}', restraint_wt=rest_wt, restraintmask=restraintmask)
    processor.add_npt_process(total_step=500000, process_name='eq_npt')
    processor.add_nvt_process(total_step=500000, process_name='eq_nvt')
    processor.add_npt_process(total_step=step_num, is_production=True, process_name='production')
    
    # Run simulation
    simulator = Simulator(processor)
    simulator.run_simulation(with_gpu)

@main.command(short_help='Automatically preparing required files for MD only.')
@click.argument('protein_file', type=click.Path(exists=True))
@click.argument('molecule_file', type=click.Path(exists=True), required=False)
@click.option('--charge', '-c', default=0, show_default=True, type=int, help='Charge of molecule.')
@click.option('--multiplicity', '-m', default=1, show_default=True, type=int, help='Multiplicity of molecule.')
@click.option('--solvent', '-s', default='water', show_default=True, type=str, help='Solvent used in molecule preparation.')
@click.option('--prefix', '-p', default=None, type=str, help='Prefix of output files.')
@click.option('--parallel', '-n', default=os.cpu_count(), show_default=True, type=int, help='Number of parallel processes.')
@click.option('--time', '-t', default=100, show_default=True, type=int, help='Total time(ns) of simulation. Default to 100 ns.')
@click.option('-bcc', is_flag=True, help='Use existing BCC charges instead of RESP.')
@click.option('--keep-cood', '-k', is_flag=True, help='Keep atoms coordinates from the input molecular file after calculating resp charge, instead of using the gaussian optimized output.')
@click.option('--keep-water', '-w', is_flag=True, help='Keep water molecules in the original protein file.')
@click.option('--box-size', '-b', default=12, show_default=True, type=float, help='TIP3P Water Box size of simulation. Default to 12 Angstrom.')
@click.option('--overwrite', '-O', is_flag=True, help='Overwrite existing files.')
def prepare(protein_file, molecule_file, charge, multiplicity, solvent, prefix, parallel, time, bcc, keep_cood, keep_water, box_size, overwrite):
    '''
    Prepare required files for MD.\n
    protein_file : Specify protein file (PDB format) path for MD.\n
    molecule_file : Specify molecule file (PDB format) path for MD. If not specified, apo system(protein only) will be simulated.\n
    '''
    from pyCADD.Dynamic import Processor
    processor = Processor() if molecule_file is not None else Processor(apo=True)
    
    if molecule_file is not None:
        if not bcc:
            processor.molecule_prepare(
                molecule_file, charge, multiplicity, parallel, solvent, overwrite=overwrite, method='resp', keep_origin_cood=keep_cood)
        else:
            processor.molecule_prepare(
                molecule_file, charge, multiplicity, parallel, solvent, overwrite=overwrite, method='bcc')
            
    processor.protein_prepare(protein_file, keep_water=keep_water)
    prefix = os.path.basename(os.getcwd()) if prefix is None else prefix
    processor.leap_prepare(prefix=prefix, box_size=box_size)
    step_num = int(time * 1000 / 0.002)
    
    water_resnum = processor.get_water_resnum()
    water_resnum_start = int(water_resnum[0])
    water_resnum_end = int(water_resnum[-1])
    
    processor.creat_minimize_input(restraint=True, restraint_mask=f"':1-{water_resnum_start-1}'", file_name="stepA.in")
    processor.creat_minimize_input(restraint=True, restraint_mask=f"':{water_resnum_start}-{water_resnum_end}'", file_name="stepB.in")
    processor.creat_minimize_input(file_name="stepC.in")
    processor.creat_heat_input(file_name="heat.in")
    # restraint md equilibration
    restraintmask = "'!(:WAT,Na+,Cl-,K+,K) & !@H= & !@H'"
    for rest_wt in [4.0, 3.5, 3.0, 2.5, 2.0, 1.0, 0]:
        processor.creat_npt_input(total_step=5000, file_name=f"eq_npt_reswt{rest_wt}.in", restraint_wt=rest_wt, restraintmask=restraintmask)
    processor.creat_npt_input(total_step=500000, file_name="eq_npt.in")
    processor.creat_nvt_input(total_step=500000, file_name="eq_nvt.in")
    processor.creat_npt_input(total_step=step_num, file_name="production.in")

@main.command(short_help='Running molecular dynamics simulation with prepared files.')
@click.argument('top_file', type=click.Path(exists=True))
@click.argument('inpcrd_file', type=click.Path(exists=True))
@click.option('--with-gpu', '-g', default=0, show_default=True, type=int, help='Specify GPU device code used in simulation. Default to 0')
def simulate(top_file, inpcrd_file, with_gpu):
    '''
    Run molecular dynamics simulation.\nOnly supports running workflows automatically generated by prepare.
    top_file : Specify solvated complex topology file (TOP format) path for MD.\n
    inpcrd_file : Specify solvated complex input coordinate file (INPCRD format) path for MD.
    '''
    from pyCADD.Dynamic import Processor, Simulator
    processor = Processor()
    processor.set_comsolvate_file(top_file, 'top')
    processor.set_comsolvate_file(inpcrd_file, 'crd')
    processor.add_process('input_file/stepA.in', 'stepA', 'minimize')
    processor.add_process('input_file/stepB.in', 'stepB', 'minimize')
    processor.add_process('input_file/stepC.in', 'stepC', 'minimize')
    processor.add_process('input_file/heat.in', 'heat')
    for rest_wt in [4.0, 3.5, 3.0, 2.5, 2.0, 1.0, 0]:
        processor.add_process(f'input_file/eq_npt_reswt{rest_wt}.in', f'eq_npt_reswt{rest_wt}', 'npt')
    processor.add_process('input_file/eq_npt.in', 'eq_npt', 'npt')
    processor.add_process('input_file/eq_nvt.in', 'eq_nvt', 'nvt')
    processor.add_process('input_file/production.in', 'production', 'npt')
    simulator = Simulator(processor)
    simulator.run_simulation(with_gpu)

@main.command(short_help='Simple post analysis for MD simulation.')
@click.option('-y', type=click.Path(exists=True), help='Trajectory file path.', prompt='Please specify trajectory file path')
@click.option('-sp', type=click.Path(exists=True), help='Solvated complex topology file path.', prompt='Please specify solvated complex topology file path')
@click.option('-lp', type=click.Path(exists=True), help='Ligand topology file path. If not specified, apo system analysis will be performed.')
@click.option('-rp', type=click.Path(exists=True), help='Receptor topology file path.', prompt='Please specify receptor topology file path')
@click.option('-cp', type=click.Path(exists=True), help='Complex topology file path.', prompt='Please specify complex topology file path')
@click.option('-ro', type=click.Path(exists=True), help='Molecular dynamics output file path.', prompt='Please specify molecular dynamics output file path')
@click.option('--no-hbond', '-nh', is_flag=True, help='Disable calculating and tracing of hydrogen bonds.')
@click.option('--no-rmsd', '-nd', is_flag=True, help='Disable calculating of RMSD.')
@click.option('--no-rmsf', '-nf', is_flag=True, help='Disable calculating of RMSF.')
# @click.option('--no-extract', '-ne', is_flag=True, help='Disable extracting of lowest energy structure.')
@click.option(
    '--decomp', '-d', 
    is_flag=True,
    help='Performing MM-GBSA energy decomposition from START_FRAME(INT1) to END_FRAME(INT2) with STEP_SIZE(INT3).')
@click.option(
    '--nmode', '-e', 
    is_flag=True,
    help='Performing entropy calculation with normal mode(nmode) from START_FRAME(INT1) to END_FRAME(INT2) with STEP_SIZE(INT3).')
@click.option('--parallel', '-n', default=os.cpu_count(), show_default=True, type=int, help='Number of parallel processes used in energy calculation.')
def analysis(y, sp, lp, rp, cp, ro, no_hbond, no_rmsd, no_rmsf, decomp, nmode, parallel=None):
    '''
    Post-process for molecular dynamics simulation.\n
    Workflow:\n
        1. Calculate RMSD\n
        2. Calculate RMSF\n
        3. Calculate and trace distance, angle and lifetime of hydrogen bonds\n
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

    if decomp:
        decomp_start_fm = input('Please specify energy decomposition START_FRAME:\n')
        decomp_end_fm = input('Please specify energy decomposition END_FRAME:\n')
        decomp_step_size = input('Please specify energy decomposition STEP_SIZE:\n')

    if nmode:
        nmode_start_fm = input('Please specify nmode START_FRAME:\n')
        nmode_end_fm = input('Please specify nmode END_FRAME:\n')
        nmode_step_size = input('Please specify nmode STEP_SIZE:\n')

    if not no_rmsd:
        analyzer.calc_rmsd()
    if not no_rmsf:
        analyzer.calc_rmsf()
    if not no_hbond:
        analyzer.calc_hbond()
    # if not no_extract:
    #     analyzer.extract_lowest_energy_st()
    if decomp:
        analyzer.creat_energy_inputfile(start_frame=decomp_start_fm, end_frame=decomp_end_fm, interval=decomp_step_size, job_type='decomp')
        analyzer.run_energy_calc(cpu_num=parallel)
        os.system('rm _MMPBSA_*')
    if nmode:
        analyzer.creat_energy_inputfile(start_frame=nmode_start_fm, end_frame=nmode_end_fm, interval=nmode_step_size, job_type='entropy')
        analyzer.run_energy_calc(cpu_num=parallel)
        os.system('rm _MMPBSA_*')

    
