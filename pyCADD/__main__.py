import os
import sys
import logging
import getopt

logger = logging.getLogger(__name__)

from pyCADD.utils.ui import UI
enter_text = '[bold]Enter the Code of Options'

from pyCADD import __version__

def main():
    if len(sys.argv) != 1:
        try:
            opts, args = getopt.getopt(sys.argv[1:], 'hvV', ['help', 'version'])
        except getopt.GetoptError as err:
            print(err)
            sys.exit(2)
        for opt, arg in opts:
            if opt in ('-v', '--version', '-V'):
                print('pyCADD Version', __version__)
                sys.exit(0)
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