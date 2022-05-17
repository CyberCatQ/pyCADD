import click

@click.command()
@click.argument('uniprot_id', type=str)
@click.option('--del_mutations', '-m', is_flag=True, default=True, help='Delete mutated(not wide type) crystal.')
@click.option('--del_ignore', '-e', is_flag=True, default=True, help='Delete small molecule which is not ligand (e.g. solvent molecules).')
@click.option('--output_format', '-o', type=click.Choice(['csv', 'in', 'ini']), default='in', help='Output format.')
def main(uniprot_id, del_mutations, del_ignore, output_format):
    '''
    Get PDB crystals info with uniprot_id.
    '''
    from pyCADD.query import QueryClient
    client = QueryClient(uniprot_id)
    client.clean_pdb_data(del_mutations, del_ignore)
    client.save(f'{uniprot_id}.{output_format}')