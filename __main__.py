import os
import sys
from pyCADD.ui import UI
enter_text = '[bold]Enter the Code of Options'

def main():
    ui = UI()
    options = [
        '1. Simple Mode',
        '2. Multiple Mode',
        '3. VSW',
        '0. Exit'
    ]
    ui.create_panel(options)
    flag = ui.get_input(enter_text, ['1', '2', '3'])

    if flag == '0':
        sys.exit(0)

    elif flag == '1':
        os.system('run python3 -m pyCADD.Dock')

    elif flag == '2':
        os.system('run python3 -m pyCADD.Multidock')

    elif flag == '3':
        os.system('run python3 -m pyCADD.VSW')

if __name__ == '__main__':
    main()