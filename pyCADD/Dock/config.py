from typing import Iterable

GLIDE_FORCEFIELD = 'OPLS3e'

class BaseConfig:
    '''
    配置项基类
    '''
    def __init__(self) -> None:
        pass

class DefaultDataConfig(BaseConfig):
    '''
    提取数据配置项
    '''
    def __init__(self, precision:str='SP') -> None:
        super().__init__()
        self.default_prop_XP = [
            'PDB', 
            'Ligand', 
            'Original', 
            'Docking_Score', 
            'MMGBSA_dG_Bind', 
            'rmsd', 
            'precision', 
            'Site_Score', 
            'Volume',
            'ligand_efficiency', 
            'XP_Hbond', 
            'rotatable_bonds', 
            'ecoul', 
            'evdw', 
            'emodel', 
            'energy', 
            'einternal', 
            'activity']
        self.default_prop_SP = [
            'PDB', 
            'Ligand', 
            'Original', 
            'Docking_Score', 
            'MMGBSA_dG_Bind', 
            'rmsd', 
            'precision', 
            'Site_Score', 
            'Volume', 
            'ligand_efficiency', 
            'rotatable_bonds', 
            'ecoul', 
            'evdw', 
            'emodel', 
            'energy', 
            'einternal',
            'lipo', 
            'hbond', 
            'metal', 
            'rewards', 
            'erotb', 
            'esite', 
            'activity'
        ]
        if precision == 'XP':
            self.properties = self.default_prop_XP
        elif precision == 'SP':
            self.properties = self.default_prop_SP
        else:
            raise ValueError('Precision must be XP or SP')

class DataConfig(DefaultDataConfig):
    def __init__(self, precision, properties:list=None) -> None:
        super().__init__(precision)
        self.properties = self.properties.extend(properties) if isinstance(properties, Iterable) else self.properties