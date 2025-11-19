import click
from pyCADD.Density.base import Gauss


@click.group()
def main():
    """
    Command line interface for pyCADD-Density.
    """
    pass


@main.command()
@click.argument("structure_file", type=click.Path(exists=True))
@click.option("--charge", "-c", default=0, show_default=True, help="Charge of the system.")
@click.option("--multiplicity", "-m", default=1, show_default=True, help="Multiplicity of the system.")
@click.option("--dft", "-d", default="B3LYP", show_default=True, help="DFT functional to use.")
@click.option("--basis_set", "-b", default="6-31G*", show_default=True, help="Basis set to use.")
@click.option("--solvent", "-s", default=None, show_default=True, help="Solvent model to use.")
@click.option("--loose", "-l", is_flag=True, help="Use loose optimization.")
@click.option("--dispersion_correct", "-D", is_flag=True, help="Use dispersion correction.")
@click.option("--parallel", "-p", type=int, default=1, show_default=True, help="Number of parallel processors.")
@click.option("--memory", "-M", type=str, default="4GB", show_default=True, help="Memory to allocate.")
@click.option("--save_dir", "-S", default=".", show_default=True, help="Directory to save results.")
def optimize(
    structure_file,
    charge,
    multiplicity,
    dft,
    basis_set,
    solvent,
    loose,
    dispersion_correct,
    parallel,
    memory,
    save_dir,
):
    """
    Optimize the molecular structure.
    """
    gauss = Gauss()
    gauss.calc_opt(
        structure_file=structure_file,
        charge=charge,
        multiplicity=multiplicity,
        dft=dft,
        basis_set=basis_set,
        solvent=solvent,
        loose=loose,
        dispersion_correct=dispersion_correct,
        cpu_num=parallel,
        mem_use=memory,
        save_dir=save_dir,
    )


@main.command()
@click.argument("structure_file", type=click.Path(exists=True))
@click.option("--charge", "-c", default=0, show_default=True, help="Charge of the system.")
@click.option("--multiplicity", "-m", default=1, show_default=True, help="Multiplicity of the system.")
@click.option("--dft", "-d", default="B3LYP", show_default=True, help="DFT functional to use.")
@click.option("--basis_set", "-b", default="6-31G*", show_default=True, help="Basis set to use.")
@click.option("--solvent", "-s", default=None, show_default=True, help="Solvent model to use.")
@click.option("--esp_calculate", "-e", is_flag=True, help="Calculate ESP.")
@click.option("--parallel", "-p", type=int, default=1, show_default=True, help="Number of parallel processors.")
@click.option("--memory", "-M", type=str, default="4GB", show_default=True, help="Memory to allocate.")
@click.option("--save_dir", "-S", default=".", show_default=True, help="Directory to save results.")
def single_point(
    structure_file,
    charge,
    multiplicity,
    dft,
    basis_set,
    solvent,
    esp_calculate,
    parallel,
    memory,
    save_dir,
):
    """
    Perform single point energy calculation.
    """
    gauss = Gauss()
    gauss.calc_energy(
        structure_file=structure_file,
        charge=charge,
        multiplicity=multiplicity,
        dft=dft,
        basis_set=basis_set,
        solvent=solvent,
        esp_calculate=esp_calculate,
        cpu_num=parallel,
        mem_use=memory,
        save_dir=save_dir,
    )


@main.command()
@click.argument("structure_file", type=click.Path(exists=True))
@click.option("--charge", "-c", default=0, show_default=True, help="Charge of the system.")
@click.option("--multiplicity", "-m", default=1, show_default=True, help="Multiplicity of the system.")
@click.option("--dft", "-d", default="B3LYP", show_default=True, help="DFT functional to use.")
@click.option("--basis_set", "-b", default="6-31G*", show_default=True, help="Basis set to use.")
@click.option("--solvent", "-s", default="water", show_default=True, help="Solvent model to use.")
@click.option("--parallel", "-p", type=int, default=1, show_default=True, help="Number of parallel processors.")
@click.option("--memory", "-M", type=str, default="4GB", show_default=True, help="Memory to allocate.")
@click.option("--save_dir", "-S", default=".", show_default=True, help="Directory to save results.")
def resp(structure_file, charge, multiplicity, dft, basis_set, solvent, parallel, memory, save_dir):
    """
    Calculate RESP charges for the molecular structure.
    """
    gauss = Gauss()
    gauss.calc_resp(
        structure_file=structure_file,
        charge=charge,
        multiplicity=multiplicity,
        dft=dft,
        basis_set=basis_set,
        solvent=solvent,
        cpu_num=parallel,
        mem_use=memory,
        save_dir=save_dir,
    )


@main.command()
@click.argument("structure_file", type=click.Path(exists=True))
@click.option("--charge", "-c", type=int, default=0, show_default=True, help="Charge of the system.")
@click.option("--multiplicity", "-m", type=int, default=1, show_default=True, help="Multiplicity of the system.")
@click.option("--dft", "-d", type=str, default="B3LYP", show_default=True, help="DFT functional to use.")
@click.option("--basis_set", "-b", type=str, default="6-31G*", show_default=True, help="Basis set to use.")
@click.option("--solvent", "-s", type=str, default="water", show_default=True, help="Solvent model to use.")
@click.option("--delta", "-D", type=float, default=0.5, show_default=True, help="Delta value for RESP2 calculation.")
@click.option("--parallel", "-p", type=int, default=1, show_default=True, help="Number of parallel processors.")
@click.option("--memory", "-M", type=str, default="4GB", show_default=True, help="Memory to allocate.")
@click.option("--save_dir", "-S", default=".", show_default=True, help="Directory to save results.")
def resp2(structure_file, charge, multiplicity, dft, basis_set, solvent, delta, parallel, memory, save_dir):
    """
    Calculate RESP2 charges for the molecular structure.
    """
    gauss = Gauss()
    gauss.calc_resp2(
        structure_file=structure_file,
        charge=charge,
        multiplicity=multiplicity,
        dft=dft,
        basis_set=basis_set,
        solvent=solvent,
        cpu_num=parallel,
        mem_use=memory,
        delta=delta,
        save_dir=save_dir,
    )
