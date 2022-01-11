import sys

from pyCADD.Dock.ui import UI_dock

if __name__ == '__main__':
    
    enter_text = '[bold]Enter the Code of Options'
    ui_dock = UI_dock()

    flag = ui_dock.get_input(enter_text, choices=[str(i) for i in range(len(ui_dock.main_options))], default='0')
    if flag == '0':
        sys.exit(0)

    ui_dock.run(flag)
