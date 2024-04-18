import os
import shutil
import unittest
from tempfile import TemporaryDirectory

from schrodinger.structure import (Structure, StructureReader, _Molecule,
                                   _Residue)
from schrodinger.structutils.analyze import Ligand

from pyCADD.Dock.schrodinger.common import (BaseMaestroFile, DockResultFile,
                                            GridFile, MaestroFile, MetaData)
from pyCADD.Dock.schrodinger.core import dock, grid_generate, minimize
from pyCADD.utils.common import ChDir

from . import DEBUG, TEST_PDB_FILE_PATH, init_logger

init_logger()


class TestMetaData(unittest.TestCase):
    def test_parse_from_filename(self):
        file_name = '1ABC_LIG_glide-dock_LIG_SP.maegz'
        metadata = MetaData.parse_from_filename(file_name)
        self.assertEqual(metadata.pdbid, '1ABC')
        self.assertEqual(metadata.ligand_name, 'LIG')
        self.assertEqual(metadata.action, 'glide-dock')
        self.assertEqual(metadata.docking_ligand_name, 'LIG')
        self.assertEqual(metadata.precision, 'SP')

    def test_parse_from_dict(self):
        metadata_dict = {
            'pdbid': '1ABC',
            'ligand_name': 'LIG',
            'action': 'glide-dock',
            'docking_ligand_name': 'LIG',
            'precision': 'SP'
        }
        metadata = MetaData.parse_from_dict(metadata_dict)
        self.assertEqual(metadata.pdbid, '1ABC')
        self.assertEqual(metadata.ligand_name, 'LIG')
        self.assertEqual(metadata.action, 'glide-dock')
        self.assertEqual(metadata.docking_ligand_name, 'LIG')
        self.assertEqual(metadata.precision, 'SP')

    def test_parse_from_metadata(self):
        metadata = MetaData(pdbid='1ABC', ligand_name='LIG',
                            action='glide-dock', docking_ligand_name='LIG', precision='SP')
        new_metadata = MetaData.parse_from_metadata(metadata)
        self.assertEqual(new_metadata.pdbid, '1ABC')
        self.assertEqual(new_metadata.ligand_name, 'LIG')
        self.assertEqual(new_metadata.action, 'glide-dock')
        self.assertEqual(new_metadata.docking_ligand_name, 'LIG')
        self.assertEqual(new_metadata.precision, 'SP')

    def test_generate_file_name(self):
        metadata = MetaData(pdbid='1ABC', ligand_name='LIG',
                            action='glide-dock', docking_ligand_name='LIG', precision='SP')
        file_name = metadata.generate_file_name()
        self.assertEqual(file_name, '1ABC_LIG_glide-dock_LIG_SP')

    def test_str_representation(self):
        metadata = MetaData(pdbid='1ABC', ligand_name='LIG',
                            action='glide-dock', docking_ligand_name='LIG', precision='SP')
        self.assertEqual(str(
            metadata), "{'pdbid': '1ABC', 'ligand_name': 'LIG', 'ligand_resnum': -1, 'action': 'glide-dock', 'docking_ligand_name': 'LIG', 'precision': 'SP'}")

    def test_repr_representation(self):
        metadata = MetaData(pdbid='1ABC', ligand_name='LIG',
                            action='glide-dock', docking_ligand_name='LIG', precision='SP')
        self.assertEqual(repr(
            metadata), "{'pdbid': '1ABC', 'ligand_name': 'LIG', 'ligand_resnum': -1, 'action': 'glide-dock', 'docking_ligand_name': 'LIG', 'precision': 'SP'}")

    def test_copy(self):
        metadata = MetaData(pdbid='1ABC', ligand_name='LIG',
                            action='glide-dock', docking_ligand_name='LIG', precision='SP')
        new_metadata = metadata.copy()
        self.assertEqual(new_metadata.pdbid, '1ABC')
        self.assertEqual(new_metadata.ligand_name, 'LIG')
        self.assertEqual(new_metadata.action, 'glide-dock')
        self.assertEqual(new_metadata.docking_ligand_name, 'LIG')
        self.assertEqual(new_metadata.precision, 'SP')
        self.assertNotEqual(id(metadata), id(new_metadata))

    def test_set_get_delete(self):
        metadata = MetaData(pdbid='1ABC', ligand_name='LIG',
                            action='glide-dock', docking_ligand_name='LIG', precision='SP')
        metadata.set('new_key', 'new_value')
        self.assertEqual(metadata.get('new_key'), 'new_value')
        self.assertIsNone(metadata.get('non_existent_key', None))
        with self.assertRaises(AttributeError):
            metadata.get('non_existent_key')
        metadata.delete('new_key')
        self.assertIsNone(metadata.get('new_key', None))


class TestBaseMaestroFile(unittest.TestCase):
    def test_init(self):
        path = TEST_PDB_FILE_PATH
        metadata = MetaData.parse_from_filename(path)
        maestro_file_metadata = BaseMaestroFile(path, metadata=metadata)
        maestro_file_filename = BaseMaestroFile(path, None)
        self.assertEqual(
            maestro_file_metadata.metadata.__dict__, metadata.__dict__)
        self.assertEqual(maestro_file_metadata.metadata.__dict__,
                         maestro_file_filename.metadata.__dict__)


class TestGridFile(unittest.TestCase):
    def test_init(self):
        path = TEST_PDB_FILE_PATH.replace('.pdb', '.zip')
        metadata = MetaData.parse_from_filename(path)
        grid_file = GridFile(path, metadata=metadata, exist=False)
        self.assertEqual(grid_file.metadata.__dict__, metadata.__dict__)


class TestMaestroFile(unittest.TestCase):
    def test_init(self):
        path = TEST_PDB_FILE_PATH
        maestro_file = MaestroFile(path)
        self.assertEqual(maestro_file.file_path, path)

    def test_st_reader(self):
        path = TEST_PDB_FILE_PATH
        maestro_file = MaestroFile(path)
        st_reader = maestro_file.st_reader
        self.assertIsInstance(st_reader, StructureReader)

    def test_structures(self):
        path = TEST_PDB_FILE_PATH
        maestro_file = MaestroFile(path)
        structures = maestro_file.structures
        self.assertIsInstance(structures, list)
        self.assertIsInstance(structures[0], Structure)
        self.assertEqual(len(structures), 1)

    def test_get_structure(self):
        path = TEST_PDB_FILE_PATH
        structure = MaestroFile.get_structure(path)
        self.assertIsInstance(structure, Structure)

    def test_get_residue(self):
        path = TEST_PDB_FILE_PATH
        residue = MaestroFile(path).get_residue(resnum=500)
        self.assertIsInstance(residue, _Residue)

    def test_get_molecule_by_res(self):
        path = TEST_PDB_FILE_PATH
        molecule = MaestroFile(path).get_molecule_by_res(resnum=500)
        self.assertIsInstance(molecule, _Molecule)

    def test_get_covalent_bond(self):
        path = TEST_PDB_FILE_PATH
        covalent_bonds = MaestroFile(path).get_covalent_bond(resnum=500)
        self.assertIsNone(covalent_bonds)

    def get_molnum_by_res(self):
        path = TEST_PDB_FILE_PATH
        molnum = MaestroFile(path).get_molnum_by_res(resnum=500)
        self.assertEqual(molnum, 4)

    def test_find_ligands(self):
        path = TEST_PDB_FILE_PATH
        ligands = MaestroFile(path).find_ligands()
        self.assertIsInstance(ligands, list)
        self.assertIsInstance(ligands[0], Ligand)
        self.assertEqual(len(ligands), 1)
        self.assertEqual(ligands[0].mol_num, 4)


class TestDockResultFile(unittest.TestCase):
    def _dock(self, include_recep: bool = False):
        testdir = os.getcwd()
        shutil.copy(TEST_PDB_FILE_PATH, testdir)
        path = os.path.join(testdir, os.path.basename(TEST_PDB_FILE_PATH))
        minimized_file = minimize(MaestroFile(path))
        lig = minimized_file.find_ligands()[0]
        lig.st.write('ligand.mae')
        grid_file = grid_generate(
            minimized_file, box_center_molnum=lig.mol_num)
        return dock(grid_file, 'ligand.mae',
                    include_receptor=include_recep,
                    save_dir=os.path.join(testdir, f"dock_include_{include_recep}"))

    def test_init(self):
        path = TEST_PDB_FILE_PATH  # Not a real dock result file
        metadata = MetaData(pdbid='3OAP', ligand_name='9cr', action='glide-dock',
                            docking_ligand_name='9CR', precision='SP')
        metadata.set('include_receptor', True)
        dock_result_file = DockResultFile(path, metadata=metadata)
        self.assertEqual(dock_result_file.file_path, path)
        self.assertTrue(dock_result_file.include_receptor)

        metadata.set('include_receptor', False)
        dock_result_file = DockResultFile(path, metadata=metadata)
        self.assertFalse(dock_result_file.include_receptor)

        dock_result_file = DockResultFile(path)
        self.assertFalse(dock_result_file.include_receptor)

        dock_result_file = DockResultFile(path, include_receptor=True)
        self.assertTrue(dock_result_file.include_receptor)

    def test_get_receptor_structure(self):
        with TemporaryDirectory() as tmpdir:
            testdir = os.path.join(os.getcwd(), 'tests') if DEBUG else tmpdir
            with ChDir(testdir):
                dock_result_file_ligand_only = self._dock(False)
                dock_result_file_with_recep = self._dock(True)
            receptor_structure = dock_result_file_ligand_only.get_receptor_structure()
            self.assertIsNone(receptor_structure)
            receptor_structure = dock_result_file_with_recep.get_receptor_structure()
            self.assertIsInstance(receptor_structure, Structure)

    def test_get_ligand_structure(self):
        with TemporaryDirectory() as tmpdir:
            testdir = os.path.join(os.getcwd(), 'tests') if DEBUG else tmpdir
            with ChDir(testdir):
                dock_result_file = self._dock()
            ligand_structure = dock_result_file.get_ligand_structure()
            self.assertIsInstance(ligand_structure, Structure)

    def test_get_raw_results(self):
        with TemporaryDirectory() as tmpdir:
            testdir = os.path.join(os.getcwd(), 'tests') if DEBUG else tmpdir
            with ChDir(testdir):
                dock_result_file = self._dock(False)
            raw_result_list = dock_result_file.get_raw_results()
            self.assertIsInstance(raw_result_list, list)
            self.assertTrue('r_i_docking_score' in raw_result_list[0])

    def test_get_results(self):
        with TemporaryDirectory() as tmpdir:
            testdir = os.path.join(os.getcwd(), 'tests') if DEBUG else tmpdir
            with ChDir(testdir):
                dock_result_file = self._dock(False)
            result_list = dock_result_file.get_results()
            self.assertIsInstance(result_list, list)
            self.assertTrue('Docking_Score' in result_list[0])
            self.assertTrue(result_list[0]['Docking_Score'])


if __name__ == '__main__':
    unittest.main()
