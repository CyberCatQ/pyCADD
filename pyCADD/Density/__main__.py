import sys
import argparse
from pyCADD.Density.ui import UI_Gauss

def main():
    args = arg_parse()
    enter_text = '[bold]Enter the Code of Options'
    if args.original_st is not None:
        original_st = args.original_st
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
        
def arg_parse():
    parser = argparse.ArgumentParser(description='pyCADD Density module for gaussian calculation')
    parser.add_argument('original_st', nargs='?', help='The initial structure file.', default=None)
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    main()