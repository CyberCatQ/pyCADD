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
            opts, args = getopt.getopt(sys.argv[1:], 'hv', ['help', 'version'])
        except getopt.GetoptError as err:
            print(err)
            sys.exit(2)
        for opt, arg in opts:
            if opt in ('-v', '--version', '-V'):
                print('pyCADD Version', __version__)
                sys.exit(0)
    ui = UI()
    options = [
        '1. Dock Mode',
        '2. VSW',
        '3. Gaussian Calculation',
        '0. Exit'
    ]
    ui.create_panel(options)
    flag = ui.get_input(enter_text, choices=[str(i) for i in range(len(options))], default='0')

    if flag == '0':
        sys.exit(0)
    elif flag in '12':
        if not ui.schrodinger_check:
            logger.error('Schrodinger platform is not installed.')
            return
        else:
            if os.system(r'run python3 -c "import pyCADD"') != 0:
                os.system('run python3 -m pip install pyCADD >/dev/null')
            
    if flag == '1':
        os.system('run python3 -m pyCADD.Dock')

    elif flag == '2':
        os.system('run python3 -m pyCADD.VSW')
    
    elif flag == '3':
        os.system('python -m pyCADD.Gauss')
        
if __name__ == '__main__':
    main()