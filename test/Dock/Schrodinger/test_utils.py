import os
import shutil
import unittest
from tempfile import TemporaryDirectory

from pyCADD.Dock.schrodinger.common import MaestroFile
from pyCADD.Dock.schrodinger.utils import (Job, collect_structures,
                                           convert_format, get_centroid,
                                           launch)
from pyCADD.utils.common import ChDir

from . import TEST_ASSETS_DIR, TEST_PDB_FILE_PATH


class TestSchrodingerUtils(unittest.TestCase):
    def test_launch(self):
        # Test launch function with a command that completes within the timeout
        with TemporaryDirectory() as tmpdir:
            shutil.copy(TEST_PDB_FILE_PATH, tmpdir)
            with ChDir(tmpdir):
                cmd = f"prepwizard -nopreprocess -noprotassign -noimpref {os.path.basename(TEST_PDB_FILE_PATH)} test.mae"
                job = launch(cmd)
                self.assertIsInstance(job, Job)
                self.assertEqual(job.Status, 'completed')
                os.remove('test.mae')

    def test_collect_structures(self):
        # Test collect_structures function
        list_file = os.path.join(TEST_ASSETS_DIR, 'structure_list.csv')
        data_dir = os.path.abspath(os.path.join(TEST_ASSETS_DIR, 'structures'))
        with TemporaryDirectory() as tmpdir:
            output_file = os.path.join(tmpdir, 'output.mae')
            collect_structures(list_file, output_file, data_dir, 'pdb')
            self.assertTrue(os.path.exists(output_file))

    def test_convert_format(self):
        # Test convert_format function with a supported format
        file_path = TEST_PDB_FILE_PATH
        to_format = 'mae'
        with TemporaryDirectory() as save_dir:
            converted_file = convert_format(file_path, to_format, save_dir)
            self.assertTrue(os.path.exists(converted_file))

    def test_get_centroid(self):
        # Test get_centroid function with a MaestroFile object
        file_path = TEST_PDB_FILE_PATH
        maestro_file = MaestroFile(file_path)
        centroid = get_centroid(maestro_file)
        self.assertIsInstance(centroid, tuple)
        self.assertEqual(len(centroid), 3)

        # Test get_centroid function with a Structure object
        structure = maestro_file.structures[0]
        centroid = get_centroid(structure)
        self.assertIsInstance(centroid, tuple)
        self.assertEqual(len(centroid), 3)

        # Test get_centroid function with a file path
        centroid = get_centroid(file_path)
        self.assertIsInstance(centroid, tuple)
        self.assertEqual(len(centroid), 3)
