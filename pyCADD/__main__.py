import os
import sys
import logging

logger = logging.getLogger(__name__)

from pyCADD.ui import UI
enter_text = '[bold]Enter the Code of Options'

def main():
    ui = UI()
    options = [
        '1. Simple Mode',
        '2. Multiple Mode',
        '3. VSW',
        '4. Gaussian Calculation',
        '5. Dance (Data Analyzer)',
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
        os.system('run python3 -m pyCADD.Multidock')

    elif flag == '3':
        os.system('run python3 -m pyCADD.VSW')
    
    elif flag == '4':
        os.system('python -m pyCADD.Gauss')
        
    elif flag == '5':
        os.system('python -m pyCADD.Dance')

if __name__ == '__main__':
    main()