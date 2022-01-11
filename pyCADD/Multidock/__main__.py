import sys

from pyCADD.Multidock.ui import UI_Multimode

if __name__ == '__main__':

    enter_text = '[bold]Enter the Code of Options'
    ui_multimode = UI_Multimode()
    while True:
        flag = ui_multimode.get_input(
            enter_text, choices=[str(i) for i in range(len(ui_multimode.main_options))], default='0')
        if flag == '0':
            sys.exit(0)
        ui_multimode.run(flag)
