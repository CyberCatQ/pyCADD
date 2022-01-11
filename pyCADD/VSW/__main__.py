import sys
import time

from pyCADD.VSW.ui import UI_VSW

if __name__ == '__main__':

    enter_text = '[bold]Enter the Code of Options'
    ui_vsw = UI_VSW()
    
    while True:
        time.sleep(0.5)
        flag = ui_vsw.get_input(enter_text, choices=[str(i) for i in range(len(ui_vsw.main_options))], default='0')
        if flag == '0':
            sys.exit(0)
        ui_vsw.run(flag)
