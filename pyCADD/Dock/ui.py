import os
from pyCADD.Dock.common import LigandFile
from pyCADD.utils.ui import UI

class UI_dock(UI):
    '''
    单晶体对接UI
    '''

    def __init__(self, menu_name: str = 'Dock Mode') -> None:
        super().__init__(menu_name=menu_name)
        self.main_options = [
            '1. PDB file download + Optimize',
            '2. PDB file download + Optimize + Generate grid file (Size 20A)', 
            '3. Generate grid file (custom Size) only', 
            '4. Redock (SP precision, Calculate MM-GBSA optionally)', 
            '5. Redock (XP precision, Calculate MM-GBSA optionally)', 
            '6. Specified ligand docking', 
            '7. ADMET Prediction of ligand',
            '0. Exit'
        ]
        self.create_panel(self.main_options)
    
    def run(self, flag):

        # 需要schrodinger安装目录中的run环境运行
        try:
            from pyCADD.Dock import Docker
        except ImportError:
            os.system('run python3 -m pip install rich ConcurrentLogHandler 1>/dev/null 2>&1')
            from pyCADD.Dock import Docker
            
        docker = Docker()
        keep_chain = self.get_confirm('Keep singel chain?')
        if flag == '1':
            if keep_chain:
                docker.keep_chain()
            docker.minimize()

        elif flag == '2' or flag == '3':
            if keep_chain:
                docker.keep_chain()
            size = 20
            docker.minimize()
            if flag == '3':
                size = self.get_input('Enter the grid size (A)', default=size)
            docker.grid_generate(gridbox_size=size)

        elif flag == '4' or flag == '5':
            precision = 'SP' if flag == '4' else 'XP'
            flag_calc_mmgbsa = self.get_confirm('Calculate MM-GBSA?')
            if keep_chain:
                docker.keep_chain()
            docker.minimize()
            docker.grid_generate()
            docker.set_precision(precision)
            docker.set_calc_rmsd(True)
            docker.split_complex()
            docker.dock()
            if flag_calc_mmgbsa:
                docker.calc_mmgbsa()

        elif flag == '6':
            precision = input('Enter the docking precision(SP|XP):').strip().upper()
            ligand_file = LigandFile(input('Enter the ligand PATH:').strip())
            flag_calc_mmgbsa = self.get_confirm('Calculate MM-GBSA?')

            docker.minimize()
            docker.grid_generate()
            docker.set_precision(precision)
            docker.set_calc_rmsd(False)
            docker.split_complex()
            docker.dock(ligand_file)
            if flag_calc_mmgbsa:
                docker.calc_mmgbsa()
        
        elif flag == '7':

            lig_file = self.get_input('Enter the ligand PATH (Default to internal ligand)').strip()
            
            lig_file = None if not lig_file else LigandFile(lig_file)
            docker.calc_admet(lig_file)
            docker.extra_admet_data()
            docker.save_admet_data()

        if flag in '456':
            docker.save_docking_data()
