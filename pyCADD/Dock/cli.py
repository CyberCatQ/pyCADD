import os

import click


@click.group()
def cli_main():
    '''
    Command line interface for pyCADD Dock Module.
    '''
    pass

@cli_main.command(short_help='Download PDB file from RCSB.')
@click.option('-i','--id', default=None, help='Single PDB ID to download.')
@click.option('-f', '--file', type=str, help='Ensemble input file containing a PDB column.')
@click.option('--save_dir', '-s', default=os.getcwd(), help='Directory to save the result file.')
@click.option('--overwrite', '-O', is_flag=True, help='Overwrite the existing file.')
def download(id:str=None, file:str=None, save_dir:str=None, overwrite:bool=False):
    from pyCADD.Dock.common import EnsembleInputFile
    from pyCADD.utils.tool import download_pdb, download_pdb_list
    if id is not None:
        download_pdb(id, save_dir, overwrite=overwrite)
    elif file is not None:
        input_file = EnsembleInputFile.parse_file(file)
        pdb_list = input_file.get_pdbid_list()
        download_pdb_list(pdb_list, save_dir=save_dir, overwrite=overwrite)


@cli_main.command(short_help='Minimize the structures before docking.')
@click.argument('file_path')
@click.option('--ph', '-p', default=7.4, type=float, help='Specify the pH value for protonation. Default 7.4.')
@click.option('--force_field', '-ff', default='OPLS4', type=click.Choice(['OPLS4', 'OPLS3e', 'OPLS3', 'OPLS_2005']), help='Force field for minimization. Default OPLS4.')
@click.option('--side_chain', '-s', is_flag=True, help='Add the missing side chain.')
@click.option('--missing_loop', '-l', is_flag=True, help='Add the missing loop.')
@click.option('--del_water', '-w', type=float, default=None, help='Delete water molecules outside the specified distance around the ligand. If specified, a float value is required.')
@click.option('--save_dir', '-s', default=os.getcwd(), help='Directory to save the result file.')
@click.option('--overwrite', '-O', is_flag=True, help='Overwrite the existing file.')
def minimize(file_path, ph, force_filed, side_chain, missing_loop, del_water, save_dir, overwrite):
    '''
    Minimize the structures before docking.
    '''
    from pyCADD.Dock.schrodinger.api import DockControl
    control = DockControl(file_path, save_path=save_dir)
    control.minimize(
        ph=ph, 
        force_field=force_filed, 
        fill_side_chain=side_chain, 
        add_missing_loop=missing_loop, 
        del_water=bool(del_water),
        watdist=del_water,
        overwrite=overwrite
        )

@cli_main.command(short_help='Generate the grid before docking.')
@click.argument('complex_file_path')
@click.option('--center', '-c', nargs=3, type=float, default=None, help='Specify the center of the grid (x, y, z).')
@click.option('--ligand', '-l', type=str, default=None, help='Residue name of the ligand defining grid center. Will be ignored if the center is specified.')
@click.option('--size', '-s', default=20, type=int, help='Grid box size (A).')
@click.option('--force_field', '-ff', default='OPLS4', type=click.Choice(['OPLS4', 'OPLS3e', 'OPLS3', 'OPLS_2005']), help='Force field for grid generation. Default OPLS4.')
@click.option('--save_dir', '-s', default=os.getcwd(), help='Directory to save the result file.')
@click.option('--overwrite', '-O', is_flag=True, help='Overwrite the existing file.')
def grid_generate(complex_file_path, center, ligand, size, force_field, save_dir, overwrite):
    '''
    Generate the grid before docking.
    '''
    from pyCADD.Dock.schrodinger.api import DockControl
    control = DockControl(save_path=save_dir)
    if center is not None:
        control.grid_generate(
            complex_file_path, 
            box_center=center, 
            force_field=force_field, 
            box_size=size, 
            overwrite=overwrite
            )
    elif ligand is not None:
        control.grid_generate(
            complex_file_path,
            box_center_resname=ligand,
            box_size=size,
            force_field=force_field,
            overwrite=overwrite
            )

@cli_main.command(short_help='Perform ligand docking.')
@click.argument('grid_file_path', type=str)
@click.argument('ligand_file_path', type=str)
@click.option('-p', '--precision', default='SP', required=False,type=click.Choice(['SP', 'XP']), help='Docking Precision (SP/XP), default SP.')
@click.option('--rmsd', '-r', is_flag=True,required=False, default=True, help='Calculate RMSD during docking.')
@click.option('--overwrite', '-O', is_flag=True, help='Overwrite the file.')
def dock(grid_file_path, ligand_file_path, precision, rmsd, overwrite):
    '''
    Perform ligand docking.  \n
    frid_file_path : Specify grid file path for ligand docking.
    ligand_file_path : Specify ligand file path for ligand docking.  
    '''
    from pyCADD.Dock.common import GridFile, LigandFile
    from pyCADD.Dock.core import dock as _dock
    grid_file = GridFile(grid_file_path)
    ligand_file = LigandFile(ligand_file_path)
    _dock(grid_file, ligand_file, precision=precision, calc_rmsd=rmsd, overwrite=overwrite)

@cli_main.command(short_help='Calculate MM-GBSA.')
@click.argument('file_path', type=str)
@click.option('--overwrite', '-O', is_flag=True, help='overwrite the file.')
def calc_mmgbsa(file_path, overwrite):
    '''
    Specify ligand file path for MM-GBSA calculation.
    '''
    from pyCADD.Dock.common import DockResultFile
    result_file = DockResultFile(file_path)
    result_file.calc_mmgbsa(overwrite=overwrite)

@cli_main.command(short_help='Calculate ADMET for ligand(s).')
@click.argument('file_path', default=None, type=str)
@click.option('--overwrite', '-O', is_flag=True, help='Overwrite the file.')
def calc_admet(file_path, overwrite):
    '''Specify ligand file path for ADMET calculation.'''
    from pyCADD.Dock.common import LigandFile
    ligand_file = LigandFile(file_path)
    ligand_file.calc_admet(overwrite=overwrite)

@cli_main.command(short_help='Ensemble Docking')
@click.argument('input_file_path', type=str)
@click.argument('library_file_path', type=str, required=False)
@click.option('--parallel', '-n', default=os.cpu_count(), type=int, help='Number of parallel processes.')
@click.option('--precision', '-p', default='SP', required=False, type=click.Choice(['SP', 'XP', 'HTVS']), help='Docking Precision (HTVS/SP/XP), default SP.')
@click.option('--del_water', '-w', is_flag=True, help='Delete all water molecules from crystal. If not specified, water molecules near the binding pocket within 5A will be kept.')
@click.option('--redock', is_flag=True, help='Redock ligands from crystals to grids or not.')
@click.option('--overwrite', '-O', is_flag=True, help='Overwrite the file.')
def ensemble_dock(input_file_path, library_file_path, parallel, precision, del_water, redock, overwrite):
    '''
    Perform ensemble docking. \n
    input_file_path : Specify input file path for ensemble docking.
    library_file_path : Specify compounds library file path for ensemble docking. If not specified, redocking ligands will be performed even without redock flag.
    '''
    from pyCADD.Dock import MultiDocker
    from pyCADD.Dock.common import LigandFile, MultiInputFile
    from pyCADD.Dock.data import (save_ensemble_docking_data,
                                  save_redocking_data)
    input_file = MultiInputFile.read_from_config(input_file_path)
    library_file = LigandFile(library_file_path) if library_file_path is not None else None
    
    console = MultiDocker(input_file)
    if library_file is not None:
        console.ligand_split(library_file, overwrite=overwrite)
    console.set_parallel_num(parallel)
    console.multi_minimize(overwrite=overwrite, del_water=del_water)
    console.multi_grid_generate(overwrite=overwrite)
    console.minimized_split()
    if library_file is not None:
        console.creat_mapping()

    if redock or library_file is None:
        console.multi_redock(precision=precision, overwrite=overwrite)
        redock_data_list = console.multi_extract_data(redock_data='only')
        save_redocking_data(redock_data_list, save_dir=console.result_save_dir)
        
    if library_file is not None:
        console.multi_dock(precision=precision, overwrite=overwrite)
        dock_data_list = console.multi_extract_data(redock_data=False, precision=precision)
        save_ensemble_docking_data(dock_data_list, save_dir=console.result_save_dir, precision=precision)
    
@cli_main.command(short_help='Extract ensemble docking data.')
@click.argument('input_file_path', type=str)
@click.argument('library_file_path', type=str)
@click.option('--parallel', '-n', default=os.cpu_count(), type=int, help='Number of parallel processes.')
@click.option('--precision', '-p', default='SP', required=False, type=click.Choice(['SP', 'XP', 'HTVS']), help='Docking Precision (HTVS/SP/XP), default SP.')
@click.option('--redock', is_flag=True, help='Redock ligands from crystals to grids or not.')
def extract_data(input_file_path, library_file_path, parallel, precision, redock):
    '''
    Extract ensemble docking data from finished docking results. \n
    input_file_path : Specify input file path for ensemble docking.\n
    library_file_path : Specify compounds library file path for ensemble docking.
    '''
    from pyCADD.Dock import MultiDocker
    from pyCADD.Dock.common import LigandFile, MultiInputFile
    from pyCADD.Dock.data import (save_ensemble_docking_data,
                                  save_redocking_data)
    input_file = MultiInputFile.read_from_config(input_file_path)
    library_file = LigandFile(library_file_path)
    console = MultiDocker(input_file)
    console.ligand_split(library_file)
    console.set_parallel_num(parallel)

    if redock:
        redock_data_list = console.multi_extract_data(redock_data='only')
        save_redocking_data(redock_data_list, save_dir=console.result_save_dir)
    dock_data_list = console.multi_extract_data(redock_data=False, precision=precision)
    save_ensemble_docking_data(dock_data_list, save_dir=console.result_save_dir, precision=precision)
    
@cli_main.command(short_help='Generate excel reports for every ligand to compare them with cocrystal ligands.')
@click.argument('input_file_path', type=str)
@click.argument('ligand_file_path', type=str)
@click.option('--parallel', '-n', default=os.cpu_count(), type=int, help='Number of parallel processes.')
@click.option('--precision', '-p', default='SP', required=False, type=click.Choice(['SP', 'XP']), help='Docking Precision (SP/XP), default SP.')
@click.option('--overwrite', '-O', is_flag=True, help='Overwrite the file.')
def generate_report(input_file_path, ligand_file_path, parallel, precision, overwrite):
    '''
    Generate excel reports for every ligand. \n
    input_file_path : Specify input file path for quick report.
    ligand_file_path : Specify ligand file path for quick report.
    '''
    from pyCADD.Dock import MultiDocker
    from pyCADD.Dock.common import LigandFile, MultiInputFile
    from pyCADD.Dock.data import Reporter, save_ensemble_docking_data
    input_file = MultiInputFile.read_from_config(input_file_path)
    ligand_file = LigandFile(ligand_file_path)
    console = MultiDocker(input_file)
    console.ligand_split(ligand_file, overwrite=overwrite)
    console.set_parallel_num(parallel)
    console.multi_minimize(overwrite=overwrite)
    console.multi_grid_generate(overwrite=overwrite)
    console.multi_redock(precision=precision, overwrite=overwrite, self_only=True)
    redock_data_list = console.multi_extract_data(redock_data='only')
    console.creat_mapping()
    console.multi_dock(precision=precision, overwrite=overwrite)
    dock_data_list = console.multi_extract_data(redock_data=False)
    save_ensemble_docking_data(dock_data_list, save_dir=console.result_save_dir)
    Reporter(input_file).generate_report(redock_data_list, dock_data_list, console.ligand_save_dir, console.result_save_dir)

@cli_main.command(short_help='Get ligand file from prediction result.')
@click.argument('predict_file', type=str)
@click.option('--output_file', '-o', type=str, default=None, help='Specify output file path.')
def extract_struc(predict_file, output_file):
    from pyCADD.Dock.common import get_predict_structure
    get_predict_structure(predict_file, output_file)