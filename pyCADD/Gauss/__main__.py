import sys

from pyCADD.Gauss.ui import UI_Gauss

if __name__ == '__main__':
    enter_text = '[bold]Enter the Code of Options'
    ui_gauss = UI_Gauss()

    while True:
        flag = ui_gauss.get_input(enter_text, choices=[str(i) for i in range(len(ui_gauss.main_options))], default='0')
        if flag == '0':
            sys.exit(0)
        
        ui_gauss.run(flag)
