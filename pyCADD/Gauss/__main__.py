import sys

from pyCADD.Gauss.ui import UI_Gauss

def main():
    enter_text = '[bold]Enter the Code of Options'
    if len(sys.argv) != 1:
        original_st = sys.argv[1]
    else:
        original_st = None
    ui_gauss = UI_Gauss(original_st=original_st)

    while True:
        flag = ui_gauss.get_input(enter_text, choices=[str(i) for i in range(len(ui_gauss.main_options))], default='0')
        if flag == '0':
            sys.exit(0)
        
        ui_gauss.run(flag)
        if flag not in '123':
            break

if __name__ == '__main__':
    main()