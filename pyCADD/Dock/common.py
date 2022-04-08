import re
import os
import logging

from schrodinger import structure as struc

logger = logging.getLogger(__name__)

def check_pdb(pdb: str):
    '''
    检查PDB ID合法性

    Return
    -------
    bool
        合法返回True 否则返回False
    '''

    if re.fullmatch(r'^\d[0-9a-zA-Z]{3,}$', pdb) is not None:
        return True
    else:
        return False

def get_input_pdbid() -> str:
    '''
    获取用户输入的PDBID并检查合法性

    Return
    ----------
    str
        PDB ID字符串
    '''

    pdbid = os.path.split(os.getcwd())[-1]

    if check_pdb(pdbid):        
        return pdbid
    else:
        while True:
            logger.info('To get PDB ID automatically, please change the name of crystal folder to PDBID.')
            pdbid = input('Input PDB ID:').strip().upper()
            logger.info('PDB ID: %s' % pdbid)
            if check_pdb(pdbid):
                return pdbid
            else:
                logger.warning('请输入正确的PDB ID!')

class BaseFile:
    '''
    基本文件类型
    '''
    def __init__(self, path) -> None:
        self.file_path = path
        self.file_name = os.path.split(path)[-1]
        self.file_dir = os.path.split(path)[0]
        self.file_ext = os.path.splitext(path)[-1]

class PDBFile(BaseFile):
    '''
    PDB文件类型
    '''
    def __init__(self, path) -> None:
        '''
        Parameter
        ----------
        path : str
            PDB文件路径
        '''
        super().__init__(path)
        self.pdbid = os.path.split(path)[-1].split('.')[0]
        
    def get_lig(self) -> list:
        '''
        从PDB文件获取配体小分子信息
        按行分割并返回列表

        Return
        ----------
        list
            配体小分子信息列表
        '''
        return os.popen("cat %s | grep -w -E ^HET | awk '{print $2}'" % self.pdbfile_path).read().splitlines()

    def get_lig_name(self) -> str:
        '''
        从PDB文件获取配体小分子名称
        按行分割并返回列表

        Return
        ----------
        str
            配体小分子名称
        '''

        lig_list = self.get_lig()

        # 去除水分子
        try:
            lig_list.remove('HOH')
        except ValueError:
            pass

        # 去除重复分子
        lig_set = set(lig_list)

        if len(lig_list) == 0:
            return None
        elif len(lig_list) == 1:
            return lig_list[0]
        else:
            logger.info('Crystal %s has multiple ligands: %s' % (self.pdbid, ' '.join(lig_set)))
            return self._input_from_list(lig_set)

    def _input_from_list(self) -> str:
        '''
        获取限制性输入的配体小分子名称(仅限于liglist内)

        Return
        ----------
        str
            配体小分子名称
        '''
        while True:
            ligname = input('Please specify ligand name:').strip().upper()
            if ligname in self.get_lig():
                return ligname
            else:
                logger.warning('Wrong ligand name, please try again.')

class MaestroFile(BaseFile):
    '''
    Maestro文件类型
    '''
    def __init__(self, path) -> None:
        super().__init__(path)
        self.structure = self.get_first_structure(path)
    
    @staticmethod
    def get_first_structure(file_path: str):
        '''
        获取Maestro文件中的第一个结构

        Parameter
        ----------
        file_path : str
            Maestro文件路径
        '''
        return next(struc.StructureReader(file_path))

    @staticmethod
    def convert_format(file_path, to_format: str) -> str:
        '''
        转换结构格式

        Parameter
        ----------
        file_path : str
            文件路径
        to_format : str
            转换后的格式
        '''
        st = MaestroFile.get_first_structure(file_path)
        prefix, suffix = os.path.splitext(file_path)
        if to_format == 'pdb':
            converted_file = prefix + '.pdb'
        elif to_format == 'mae':
            converted_file = prefix + '.mae'
        elif to_format == 'maegz':
            converted_file = prefix + '.maegz'
        else:
            raise ValueError('Unsupported format: %s' % to_format)

        st.write(converted_file)
        return converted_file

    def _get_mol_obj(self, ligname: str):   
        '''
        获取结构中的配体所在Molecule object

        Parameter
        ----------
        ligname : str
            配体名称
        '''
        residues = self.structure.residue
        for res in residues:
            if res.pdbres.strip() == '%s' % ligname:
                molnum = res.molecule_number
                yield self.structure.molecule[molnum]

    def _del_covalent_bond(self, ligname):
        '''
        删除共价键

        Parameter
        ----------
        ligname : str
            配体名称
        '''
        bonds = self.structure.bond

        for bond in bonds:
            resname1 = bond.atom1.getResidue().pdbres.strip()
            resname2 = bond.atom2.getResidue().pdbres.strip()
            if resname1 == '%s' % ligname or resname2 == '%s' % ligname:
                if resname1 != resname2:
                    bond_to_del = bond
                
        if not bond_to_del:
            raise RuntimeError('Can not delete covalent bonds automatically.')
                
        self.structure.deleteBond(bond_to_del.atom1, bond_to_del.atom2)
        
    def get_lig_molnum(self, ligname:str) -> str:
        '''
        以ligname为KEY 查找Maestro文件中的Molecule Number
        
        Parameter
        ----------
        ligname : str
            配体小分子名称

        Return
        ----------

        '''
        mol = next(self._get_mol_obj(ligname))

        # 判断该molecule是否仅包括小分子本身(是否存在共价连接) 自动移除共价连接
        if len(mol.residue) != 1:  
            logger.info('%s in %s : A covalent bond may exist between the ligand and residue.' % (ligname, self.file_name))
            logger.info('An attempt will be made to remove the covalent bond automatically.')
            self._del_covalent_bond(ligname)
            mol = next(self._get_mol_obj(ligname))

        return mol.number


class ReceptorInputFile(BaseFile):
    '''
    pyCADD受体列表输入文件
    '''
    def __init__(self, path) -> None:
        '''
        Parameter
        ----------
        path : str
            受体列表输入文件路径
        '''
        super().__init__(path)
        self.pdbid_list = None
        self.pairs_list = None


    def parse_file(self, file_path: str) -> None:
        '''
        解析受体输入文件
        文件中含有多行的 逗号分隔的PDBID与配体ID
        
            example:
                3A9E,REA
                3KMZ,EQO
                3KMR,EQN
                4DQM,LUF
                1DKF,BMS
                5K13,6Q7

        Parameter
        ----------
        file_path : str
            受体列表文件路径
        '''
        with open(file_path, 'r') as f:
            raw_list = f.read().splitlines()
        
        self.pairs_list = [(pdbid, ligid) for pdbid, ligid in [line.split(',') for line in raw_list]]
        self.pdbid_list = [pdbid for pdbid, ligid in self.pairs_list]
    
    def get_pairs_list(self) -> list:
        '''
        获取受体列表
        '''
        return self.pairs_list
    
    def get_pdbid_list(self) -> list:
        '''
        获取受体列表中的PDBID列表
        '''
        return self.pdbid_list
    
    def get_pdbfile_path_list(self, pdb_dir: str) -> list:
        '''
        获取受体列表中的PDB文件路径列表

        Parameter
        ----------
        pdb_dir : str
            PDB文件所在目录
        '''
        return [os.path.join(pdb_dir, pdbid + '.pdb') for pdbid in self.pdbid_list]
        
    def get_gridfile_path_list(self, grid_dir: str) -> list:
        '''
        获取受体列表中的Grid文件路径列表

        Parameter
        ----------
        grid_dir : str
            Grid文件所在目录
        '''
        return [os.path.join(grid_dir, '%s_glide_grid_%s.zip' % (pdbid, ligid)) for pdbid, ligid in self.pairs_list]