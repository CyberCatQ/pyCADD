from pyCADD.Dynamic import Processor, Simulator
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
    processor = Processor()
    processor.protein_prepare(protein_file)
    processor.molecule_prepare(
        molecule_file, charge, multiplicity, parallel, solvent)
    prefix = os.path.basename(os.getcwd()) if prefix is None else prefix
    processor.leap_prepare(prefix)
    simulator = Simulator(processor)
    simulator.creat_input_file()
    simulator.run_simulation(with_gpu)
