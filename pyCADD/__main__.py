import os
import sys
import logging

logger = logging.getLogger(__name__)

from pyCADD.utils.ui import UI
enter_text = '[bold]Enter the Code of Options'

from pyCADD import __version__

def main():
    import argparse
    parser = argparse.ArgumentParser(
        description='pyCADD: A Python-based Computer-Aided Drug Design Toolkit.',
        usage='''
    pycadd-dock: Docking module of pyCADD.
    pycadd-density: Gaussian calculation module of pyCADD.
    pycadd-dynamic: Molecular dynamics simulation module of pyCADD.
    pycadd-demand: Query module for the PDB data from Uniprot ID.

    Use CLI_COMMAND + -h or --help for more information. ''')
    parser.add_argument('-v', '-V', '--version', action='version', version='pyCADD Version ' + __version__)
    args = parser.parse_args()
    
    ui = UI()
    options = [
        '1. Query the PDB data from Uniprot ID',
        '2. Dock Mode',
        '3. VSW',
        '4. Gaussian Calculation',
        '0. Exit'
    ]
    ui.create_panel(options)
    flag = ui.get_input(enter_text, choices=[str(i) for i in range(len(options))], default='0')

    if flag == '0':
        sys.exit(0)
    elif flag in '23':
        if not ui.schrodinger_check:
            logger.error('Schrodinger platform is not installed.')
            return
        else:
            if os.system(r'run python3 -c "import pyCADD"') != 0:
                os.system('run python3 -m pip install pyCADD >/dev/null')
            
    if flag == '2':
        os.system('run python3 -m pyCADD.Dock')

    elif flag == '3':
        os.system('run python3 -m pyCADD.VSW')
    
    elif flag == '4':
        os.system('python -m pyCADD.Density')
    
    elif flag == '1':
        uniprot_id = input('Enter the Uniprot ID: ')
        os.system(f'python -m pyCADD.Demand {uniprot_id}')

if __name__ == '__main__':
    main()