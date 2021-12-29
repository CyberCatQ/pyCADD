import logging

from pyCADD.ui import UI

try:
    from pyCADD.Multidock.base import Multidock
except ImportError:
    import os
    os.system('run python3 -m pip install rich ConcurrentLogHandler')
    from pyCADD.Multidock.base import Multidock

logger = logging.getLogger('pyCADD.Multidock.UI')


class UI_Multimode(UI):
    '''
    多晶体对接UI
    '''

    def __init__(self, menu_name: str = 'Multiple Mode') -> None:
        super().__init__(menu_name=menu_name)
        self.receptor = ''
        self.ligand = ''
        self.multidocker = Multidock()

    def run(self, flag):
        # 需要schrodinger安装目录中的run环境运行
        
        if flag == '1':
            self.receptor = input('Enter the receptor list file PATH: ')
            self.multidocker.read_receptor(self.receptor)
            self.create_panel(
                additional_info='[bright_cyan]Loaded receptor file: %s' % self.receptor)

        elif flag == '2':
            self.ligand = input('Enter the ligand structure file PATH: ')
            self.multidocker.read_ligands(self.ligand)
            self.create_panel(
                additional_info='[bright_cyan]Loaded ligand file: %s' % self.ligand)

        if flag in '345678':
            if not self.receptor:
                self.create_panel()
                logger.error('No receptor loaded.')
                return
        if flag in '678':
            if not self.ligand:
                self.create_panel()
                logger.error('No ligand loaded.')
                return

        if flag == '3':
            self.multidocker.optimize()
            self.multidocker.split()
            self.create_panel()

        elif flag == '4':
            self.multidocker.grid_generate()
            self.create_panel()

        elif flag == '5':
            self.multidocker.self_dock()
            self.create_panel()

        elif flag == '6':
            self.multidocker.optimize()
            self.multidocker.split()
            self.multidocker.grid_generate()
            self.multidocker.self_dock()
            self.multidocker.map()
            self.multidocker.ensemble_dock()
            self.multidocker.save_data()
            self.create_panel()

        elif flag == '7':
            self.multidocker.optimize()
            self.multidocker.split()
            self.multidocker.grid_generate()
            self.multidocker.self_dock()
            self.multidocker.map()
            self.multidocker.ensemble_dock('XP')
            self.multidocker.save_data('XP')
            self.create_panel()

        elif flag == '8':
            precision = self.get_input('Enter the docking precision:', ['SP', 'XP'], 'SP')
            self.multidocker.map()

            additional_col = []
            if self.get_confirm('Additional columns need to be extracted?', default=False):
                additional_col = input('Enter the ORIGINAL name of columns (separated by commas): ').strip().split(',')

            self.multidocker.save_data(precision, additional_col)
            self.create_panel()


if __name__ == '__main__':

    enter_text = '[bold]Enter the Code of Options'
    ui_multimode = UI_Multimode()
    options = [
        '1. Read receptor list from file (txt/csv)',
        '2. Read structure of ligands from file (mae/maegz)',
        '3. Optimize structures of receptors',
        '4. Create grid files for receptors',
        '5. Re-dock receptors with self-ligand',
        '6. Run standard workflow of ensemble docking (SP precision)',
        '7. Run standard workflow of ensemble docking (XP precision)',
        '8. Extract and save data',
        '0. Exit'
    ]

    ui_multimode.create_panel(options)

    while True:
        flag = ui_multimode.get_input(enter_text, choices=[str(i) for i in range(len(options))], default='0')
        if flag == '0':
            break
        ui_multimode.run(flag)
