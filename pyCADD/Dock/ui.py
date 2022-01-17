from pyCADD.ui import UI

class UI_dock(UI):
    '''
    单晶体对接UI
    '''

    def __init__(self, menu_name: str = 'Simple Mode') -> None:
        super().__init__(menu_name=menu_name)
        self.main_options = [
            '1. PDB file download + Optimization',
            '2. PDB file download + Optimization + Generate grid file (Size 20A)', 
            '3. Generate grid file (custom Size) only', 
            '4. Internal ligand docking automatically (SP precision, Calculate MM-GBSA)', 
            '5. Internal ligand docking automatically (XP precision, Calculate MM-GBSA)', 
            '6. Specified ligand docking', 
            '7. Internal ligand docking automatically (SP precision)', 
            '8. Internal ligand docking automatically (XP precision)',
            '9. ADMET Prediction of ligand',
            '0. Exit'
        ]
        self.create_panel(self.main_options)
    
    def run(self, flag):

        # 需要schrodinger安装目录中的run环境运行
        try:
            from pyCADD.Dock.base import Docker
        except ImportError:
            import os
            os.system('run python3 -m pip install rich ConcurrentLogHandler')
            from pyCADD.Dock.base import Docker
            
        docker = Docker()
        if flag == '1':
            docker.minimize()

        elif flag == '2':
            docker.minimize()
            docker.grid_generate()

        elif flag == '3':
            size = input('Enter the grid box size(A):')
            docker.grid_generate(size)

        elif flag == '4':
            docker.minimize()
            docker.grid_generate()
            docker.set_precision('SP')
            docker.set_calc_rmsd(True)
            docker.split()
            docker.dock()
            docker.cal_volume()
            docker.cal_mmgbsa()

        elif flag == '5':
            docker.minimize()
            docker.grid_generate()
            docker.set_precision('XP')
            docker.set_calc_rmsd(True)
            docker.split()
            docker.dock()
            docker.cal_volume()
            docker.cal_mmgbsa()

        elif flag == '6':
            
            precision = input('Enter the docking precision(SP|XP):').strip().upper()
            ligand_file = input('Enter the ligand PATH:').strip()
            docker.minimize()
            docker.split()
            docker.grid_generate()
            docker.set_precision(precision)
            docker.dock(ligand_file)
            docker.cal_volume()
            docker.cal_mmgbsa()
            
        elif flag == '7':
            docker.minimize()
            docker.grid_generate()
            docker.set_precision('SP')
            docker.set_calc_rmsd(True)
            docker.split()
            docker.dock()
            docker.cal_volume()

        elif flag == '8':
            docker.minimize()
            docker.grid_generate()
            docker.set_precision('XP')
            docker.set_calc_rmsd(True)
            docker.split()
            docker.dock()
            docker.cal_volume()

        elif flag == '9':
            lig_file = input('Enter the ligand PATH (Default to internal ligand):').strip()
            if not lig_file:
                lig_file = None
            docker.cal_admet(lig_file)
            docker.extra_admet_data()
            docker.save_admet_data()

        if flag in '45678':
            docker.extra_data()
            docker.save_data()
