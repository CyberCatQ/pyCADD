import os
import shutil
import unittest
from tempfile import TemporaryDirectory, mkdtemp

from pyCADD.Dock.schrodinger.api import DockControl, DockEnsemble
from pyCADD.Dock.schrodinger.common import MaestroFile

from . import DEBUG, TEST_ASSETS_DIR, TEST_PDB_FILE_PATH


class TestDockControl(unittest.TestCase):

    def test_download_pdb(self):
        with TemporaryDirectory() as tmp_dir:
            testdir = os.path.join(
                os.getcwd(), 'tests') if DEBUG else tmp_dir
            control = DockControl(
                pdbid=os.path.basename(TEST_PDB_FILE_PATH).split('.')[0],
                save_path=testdir
                )

    def test_minimize(self):
        with TemporaryDirectory() as tmp_dir:
            testdir = os.path.join(
                os.getcwd(), 'tests') if DEBUG else tmp_dir
            control = DockControl(protein_file=TEST_PDB_FILE_PATH, save_path=testdir)
            minimized_file = control.minimize()
            self.assertTrue(os.path.exists(minimized_file.file_path))
    
    def test_grid_generate(self):
        with TemporaryDirectory() as tmp_dir:
            testdir = os.path.join(
                os.getcwd(), 'tests') if DEBUG else tmp_dir
            control = DockControl(protein_file=TEST_PDB_FILE_PATH, save_path=testdir)
            minimized_file = control.minimize()
            grid_file = control.grid_generate(minimized_file, box_center_resname='9CR')
            self.assertTrue(os.path.exists(grid_file.file_path))

    def test_local_workflow(self):
        with TemporaryDirectory() as tmp_dir:
            testdir = os.path.join(
                os.getcwd(), 'tests') if DEBUG else tmp_dir
            control = DockControl(protein_file=TEST_PDB_FILE_PATH, save_path=testdir)
            control.minimize()
            control.grid_generate()
            pro_file, lig_file = control.split()
            control.dock(lig_file)

class TestDockEnsemble(unittest.TestCase):

    def setUp(self):
        self.save_path = os.path.join(os.getcwd(), 'tests', 'dock_ensemble') if DEBUG else mkdtemp()
        self.test_input_file = os.path.join(TEST_ASSETS_DIR, 'ensemble', 'P10275.yml')
        self.library_file = MaestroFile(os.path.join(TEST_ASSETS_DIR, 'structures', 'multiple_ligands.sdf'))
        self.control = DockEnsemble(self.test_input_file, save_path=self.save_path, cpu_num=6)
    
    def test_load_library(self):
        self.control.load_library(self.library_file)
        for i in range(1, 4):
            self.assertTrue(
                os.path.exists(
                    os.path.join(self.save_path, "ligands", f"{self.library_file.file_prefix}-{i}.maegz")
                    )
                )
    
    def test_workflow(self):
        for pdbid in self.control.pdbid_list:
            self.assertTrue(os.path.exists(os.path.join(self.save_path, 'pdb', f"{pdbid}.pdb")))
        self.control.load_library(self.library_file)
        self.control.keep_single_chain()
        self.control.minimize()
        self.control.grid_generate()
        self.control.dock(retrospective=True)
    
    def tearDown(self):
        if not DEBUG:
            shutil.rmtree(self.save_path)