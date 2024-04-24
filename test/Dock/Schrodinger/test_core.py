import os
import shutil
import unittest
from tempfile import TemporaryDirectory

from schrodinger.structure import Structure

from pyCADD.Dock.schrodinger.common import (DockResultFile, GridFile,
                                            MaestroFile)
from pyCADD.Dock.schrodinger.core import (dock, grid_generate, minimize,
                                          split_complex)
from pyCADD.utils.common import ChDir

from . import DEBUG, TEST_PDB_FILE_PATH, init_logger

init_logger()


class TestSchrodingerCore(unittest.TestCase):
    def setUp(self):
        self.pdbid = os.path.basename(TEST_PDB_FILE_PATH).split('.')[0]
        
    def test_split_complex(self):
        # Test case 1: Splitting complex with ligand_atom_indexes
        path = TEST_PDB_FILE_PATH
        maestro_file = MaestroFile(path)
        ligand = maestro_file.find_ligands()[0]
        receptor, ligand = split_complex(
            maestro_file, ligand_atom_indexes=ligand.atom_indexes)
        self.assertIsInstance(receptor, Structure)
        self.assertIsInstance(ligand, Structure)

        # Test case 2: Splitting complex with ligand_asl
        receptor, ligand = split_complex(
            MaestroFile(path), ligand_asl='res.ptype "9CR "')
        self.assertIsInstance(receptor, Structure)
        self.assertIsInstance(ligand, Structure)

        # Test case 3: Splitting complex without providing ligand_atom_indexes or ligand_asl
        with self.assertRaises(ValueError):
            split_complex(MaestroFile(path))
        # Test case 4: Splitting complex with ligand_asl identifying multiple ligands
        with self.assertRaises(ValueError):
            split_complex(MaestroFile(path), ligand_asl='res.ptype "LYS "')

    def test_minimize_grid_dock(self):
        with TemporaryDirectory() as tmp_dir:
            testdir = os.path.join(
                os.getcwd(), 'tests') if DEBUG else tmp_dir
            os.makedirs(testdir, exist_ok=True)

            shutil.copy(TEST_PDB_FILE_PATH, testdir)
            path = os.path.join(testdir, os.path.basename(TEST_PDB_FILE_PATH))
            with ChDir(testdir):
                
                minimized_file = minimize(MaestroFile(path))
                self.assertEqual(minimized_file.file_prefix, f"{self.pdbid}__minimized")
                self.assertIsInstance(minimized_file, MaestroFile)
                
                lig = minimized_file.find_ligands()[0]
                lig.st.write('ligand.mae')
                
                with self.assertRaises(ValueError):
                    grid_generate(minimized_file, overwrite=True)
                grid_file = grid_generate(
                    minimized_file, box_center_molnum=lig.mol_num)
                self.assertIsInstance(grid_file, GridFile)
                self.assertEqual(grid_file.file_prefix, f"{self.pdbid}__glide-grid")
                
                dock_result = dock(grid_file, 'ligand.mae')
                self.assertIsInstance(dock_result, DockResultFile)
                self.assertTrue(
                    'r_i_docking_score' in dock_result.get_raw_results()[0])
                self.assertEqual(dock_result.file_prefix, f"{self.pdbid}__glide-dock__SP")


if __name__ == '__main__':
    unittest.main()
