import click

@click.command()
@click.argument('uniprot_id', type=str)
@click.option('--not_del_mutations', '-m', is_flag=True, help='DONT Delete mutated(not wide type) crystal.')
@click.option('--not_del_ignore', '-e', is_flag=True, help='DONT Delete small molecule which is not ligand (e.g. solvent molecules).')
@click.option('--cutoff', '-c', type=float, default=None, help='Cutoff of resolution.')
@click.option('--output_format', '-o', type=click.Choice(['csv', 'in', 'ini', 'yml', 'yaml']), default='yml', help='Output format.')
def main(uniprot_id, not_del_mutations, not_del_ignore, cutoff, output_format):
    '''
    Get PDB crystals info with uniprot_id.
    '''
    from pyCADD.query import QueryClient
    client = QueryClient(uniprot_id)
    del_mutations = not not_del_mutations
    del_ignore = not not_del_ignore
    client.clean_pdb_data(del_mutations, del_ignore, cutoff)
    client.save(f'{uniprot_id}.{output_format}')