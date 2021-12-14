import os
from pyCADD.ui import UI
enter_text = '[bold]Enter the Code of Options'

def main():
    ui = UI()
    options = [
        '1. Simple Mode',
        '2. Multiple Mode'
    ]
    ui.create_panel(options)
    flag = ui.get_input(enter_text, ['1', '2'])

    if flag == '1':
        os.system('run python3 -m pyCADD.Dock')

    elif flag == '2':
        os.system('run python3 -m pyCADD.Multidock')

if __name__ == '__main__':
    main()