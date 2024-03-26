import click

@click.command()
@click.argument('uniprot_id', type=str, required=False)
@click.option('--cutoff', '-c', type=float, default=None, help='Cutoff of resolution.')
@click.option('--output_format', '-o', type=click.Choice(['csv', 'in', 'ini', 'yml', 'yaml']), default='yml', help='Output format.')
@click.option('-f', '--pdb_list_file', type=click.Path(exists=True), help='CSV file with pdb list to query.')
@click.option('--pdb_column', '-p', type=str, default='PDB', help="Column name of pdb id in pdb list file. Default by 'PDBID' if not specified.")
@click.option('--generate', '-g', is_flag=True, help='Generate input file for pyCADD-Dock.')
@click.option('--not_del_mutations', '-m', is_flag=True, help='DONT Delete mutated(not wide type) crystal when generate input file for pyCADD-Dock.')
@click.option('--not_del_ignore', '-e', is_flag=True, help='DONT Delete small molecule which is not ligand (e.g. solvent molecules) when generate input file for pyCADD-Dock')
def main(uniprot_id, not_del_mutations, not_del_ignore, cutoff, output_format, pdb_list_file, pdb_column, generate):
    '''
    Get PDB crystals info with an uniprot_id, or a csv file including specified pdb id list.
    '''
    from pyCADD.Demand import QueryClient
    client = QueryClient(uniprot_id=uniprot_id, pdb_list_file=pdb_list_file)
    client.query()
    if generate:
        del_mutations = not not_del_mutations
        del_ignore = not not_del_ignore
        client.clean_pdb_data(del_mutations, del_ignore, cutoff)
        client.generate_inputfile(f'{uniprot_id}.{output_format}')