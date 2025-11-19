import os
from pyCADD.utils.ui import UI


class UI_dock(UI):

    def __init__(self, menu_name: str = "Dock Mode") -> None:
        super().__init__(menu_name=menu_name)
        self.main_options = [
            "1. PDB file download + Optimize",
            "2. PDB file download + Optimize + Generate grid file (Size 20A)",
            "3. Generate grid file (custom Size) only",
            "4. Redock (SP precision, Calculate MM-GBSA optionally)",
            "5. Redock (XP precision, Calculate MM-GBSA optionally)",
            "6. Specified ligand docking",
            "7. ADMET Prediction of ligand",
            "0. Exit",
        ]
        self.create_panel(self.main_options)

    def run(self, flag):
        raise NotImplementedError("Dock Module UI is not implemented yet.")
