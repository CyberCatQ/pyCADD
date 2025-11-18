import os

import click


@click.group()
def cli_main():
    """
    Command line interface for pyCADD Dock Module.
    """
    pass


@cli_main.command(short_help="Download PDB file from RCSB.")
@click.option("-i", "--id", default=None, help="Single PDB ID to download.")
@click.option("-f", "--file", type=str, help="Ensemble input file containing a PDB column.")
@click.option("--save_dir", "-s", default=os.getcwd(), help="Directory to save the result file.")
@click.option("--overwrite", "-O", is_flag=True, help="Overwrite the existing file.")
def download(id: str = None, file: str = None, save_dir: str = None, overwrite: bool = False):
    from pyCADD.Dock.common import EnsembleInputFile
    from pyCADD.utils.tool import download_pdb, download_pdb_list

    if id is not None:
        download_pdb(id, save_dir, overwrite=overwrite)
    elif file is not None:
        input_file = EnsembleInputFile.parse_file(file)
        pdb_list = input_file.get_pdbid_list()
        download_pdb_list(pdb_list, save_dir=save_dir, overwrite=overwrite)


@cli_main.command(short_help="Minimize the structures before docking.")
@click.argument("file_path")
@click.option(
    "--ph", "-p", default=7.4, type=float, help="Specify the pH value for protonation. Default 7.4."
)
@click.option(
    "--force_field",
    "-f",
    default="OPLS4",
    type=click.Choice(["OPLS4", "OPLS3e", "OPLS3", "OPLS_2005"]),
    help="Force field for minimization. Default OPLS4.",
)
@click.option("--side_chain", "-s", is_flag=True, help="Add the missing side chain.")
@click.option("--missing_loop", "-l", is_flag=True, help="Add the missing loop.")
@click.option(
    "--del_water",
    "-w",
    type=float,
    default=None,
    help="Delete water molecules outside the specified distance around the ligand. If specified, a float value is required.",
)
@click.option("--save_dir", "-s", default=os.getcwd(), help="Directory to save the result file.")
@click.option("--overwrite", "-O", is_flag=True, help="Overwrite the existing file.")
def minimize(file_path, ph, force_filed, side_chain, missing_loop, del_water, save_dir, overwrite):
    """
    Minimize the structures before docking.
    """
    from pyCADD.Dock.schrodinger.api import DockControl

    control = DockControl(save_path=save_dir)
    control.minimize(
        structure_file=file_path,
        ph=ph,
        force_field=force_filed,
        fill_side_chain=side_chain,
        add_missing_loop=missing_loop,
        del_water=bool(del_water),
        watdist=del_water,
        overwrite=overwrite,
    )


@cli_main.command(short_help="Generate the grid before docking.")
@click.argument("complex_file_path")
@click.option(
    "--center",
    "-c",
    nargs=3,
    type=float,
    default=None,
    help="Specify the center of the grid (x, y, z).",
)
@click.option(
    "--ligand",
    "-l",
    type=str,
    default=None,
    help="Residue name of the ligand defining grid center. Will be ignored if the center is specified.",
)
@click.option("--size", "-s", default=20, type=int, help="Grid box size (A).")
@click.option(
    "--force_field",
    "-f",
    default="OPLS4",
    type=click.Choice(["OPLS4", "OPLS3e", "OPLS3", "OPLS_2005"]),
    help="Force field for grid generation. Default OPLS4.",
)
@click.option("--save_dir", "-s", default=os.getcwd(), help="Directory to save the result file.")
@click.option("--overwrite", "-O", is_flag=True, help="Overwrite the existing file.")
def grid_generate(complex_file_path, center, ligand, size, force_field, save_dir, overwrite):
    """
    Generate the grid before docking.
    """
    from pyCADD.Dock.schrodinger.api import DockControl

    control = DockControl(save_path=save_dir)
    if center is not None:
        control.grid_generate(
            structure_file=complex_file_path,
            box_center=center,
            force_field=force_field,
            box_size=size,
            overwrite=overwrite,
        )
    elif ligand is not None:
        control.grid_generate(
            structure_file=complex_file_path,
            box_center_resname=ligand,
            box_size=size,
            force_field=force_field,
            overwrite=overwrite,
        )


@cli_main.command(short_help="Perform ligand docking.")
@click.argument("grid_file_path", type=str)
@click.argument("ligand_file_path", type=str)
@click.option(
    "-p",
    "--precision",
    default="SP",
    required=False,
    type=click.Choice(["SP", "XP", "HTVS"]),
    help="Docking Precision (SP/XP/HTVS), default SP.",
)
@click.option(
    "--force_field",
    "-f",
    default="OPLS4",
    type=click.Choice(["OPLS4", "OPLS3e", "OPLS3", "OPLS_2005"]),
    help="Force field for grid generation. Default OPLS4.",
)
@click.option("--rmsd", "-r", is_flag=True, default=True, help="Calculate RMSD during docking.")
@click.option("--overwrite", "-O", is_flag=True, help="Overwrite the file.")
@click.option("--save_dir", "-s", default=".", help="Directory to save the result file.")
def dock(grid_file_path, ligand_file_path, precision, force_field, rmsd, save_dir, overwrite):
    """
    Perform ligand docking.  \n
    frid_file_path : Specify grid file path for ligand docking.
    ligand_file_path : Specify ligand file path for ligand docking.
    """
    from pyCADD.Dock.schrodinger.api import DockControl

    control = DockControl(save_path=save_dir)
    control.dock(
        ligand_file=ligand_file_path,
        grid_file=grid_file_path,
        precision=precision,
        force_field=force_field,
        calc_rmsd=rmsd,
        overwrite=overwrite,
    )


@cli_main.command(short_help="Ensemble Docking")
@click.argument("input_file_path", type=str)
@click.argument("library_file_path", type=str, required=False)
@click.option(
    "--parallel", "-n", default=os.cpu_count(), type=int, help="Number of parallel processes."
)
@click.option(
    "--del_water",
    "-w",
    is_flag=True,
    help="Delete all water molecules from crystal. If not specified, water molecules near the binding pocket within 5A will be kept.",
)
@click.option(
    "-p",
    "--precision",
    default="SP",
    required=False,
    type=click.Choice(["SP", "XP", "HTVS"]),
    help="Docking Precision (SP/XP), default SP.",
)
@click.option(
    "--force_field",
    "-f",
    default="OPLS4",
    type=click.Choice(["OPLS4", "OPLS3e", "OPLS3", "OPLS_2005"]),
    help="Force field for grid generation. Default OPLS4.",
)
@click.option("--redock", is_flag=True, help="Redock ligands from crystals to grids or not.")
@click.option("--rmsd", "-r", is_flag=True, default=True, help="Calculate RMSD during docking.")
@click.option("--save_dir", "-s", default=".", help="Directory to save the result file.")
@click.option("--overwrite", "-O", is_flag=True, help="Overwrite the file.")
def ensemble_dock(
    input_file_path,
    library_file_path,
    parallel,
    precision,
    del_water,
    force_field,
    redock,
    rmsd,
    save_dir,
    overwrite,
):
    """
    Perform ensemble docking. \n
    input_file_path : Specify input file path for ensemble docking.
    library_file_path : Specify compounds library file path for ensemble docking. If not specified, redocking ligands will be performed even without redock flag.
    """
    from pyCADD.Dock.schrodinger.api import DockEnsemble

    control = DockEnsemble(
        input_file=input_file_path,
        save_path=save_dir,
        cpu_num=parallel,
    )
    if library_file_path is not None:
        control.load_library(library_file_path, overwrite=overwrite)
    control.keep_single_chain(overwrite=overwrite)
    control.minimize(force_field=force_field, overwrite=overwrite, del_water=del_water)
    control.grid_generate(overwrite=overwrite)
    control.dock(
        retrospective=bool(redock or not library_file_path),
        force_field=force_field,
        precision=precision,
        calc_rmsd=rmsd,
        overwrite=overwrite,
    )
