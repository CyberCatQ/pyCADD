import sys

from pyCADD.Dance.ui import UI_Dance

enter_text = '[bold]Enter the Code of Options'

if __name__ == '__main__':
    ui_dance = UI_Dance()
    ui_dance.create_panel(ui_dance.main_options)

    while True:
        flag = ui_dance.get_input(enter_text, choices=[str(
            i) for i in range(len(ui_dance.main_options))], default='0')
        if flag == '0':
            sys.exit(0)

        ui_dance.run(flag)
