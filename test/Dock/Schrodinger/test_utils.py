import os
import shutil
import unittest
from tempfile import TemporaryDirectory
from unittest.mock import patch

from pyCADD.Dock.schrodinger.common import MaestroFile
from pyCADD.Dock.schrodinger.utils import (Job, collect_structures,
                                           convert_format, get_centroid,
                                           launch)
from pyCADD.utils.common import ChDir

from . import TEST_ASSETS_DIR, TEST_PDB_FILE_PATH, init_logger

init_logger()


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

    def test_convert_format_existing_file(self):
        # Test convert_format with an existing input file
        with TemporaryDirectory() as tmpdir:
            input_file_path = TEST_PDB_FILE_PATH
            output_file_path = os.path.join(tmpdir, 'output_file.mae')
            converted_file_path = convert_format(
                input_file_path, output_file_path, save_dir=tmpdir)

            self.assertEqual(converted_file_path, os.path.join(
                tmpdir, 'output_file.mae'))
            self.assertTrue(os.path.exists(converted_file_path))

    def test_convert_format_non_existing_file(self):
        # Test convert_format with a non-existing input file
        with TemporaryDirectory() as tmpdir:
            input_file_path = 'non_existing_file.pdb'
            output_file_path = os.path.join(tmpdir, 'output_file.mae')
            save_dir = tmpdir
            with self.assertRaises(FileNotFoundError):
                convert_format(input_file_path, output_file_path, save_dir=save_dir)

    def test_convert_format_existing_output_file(self):
        # Test convert_format with an existing output file
        with TemporaryDirectory() as tmpdir:
            input_file_path = TEST_PDB_FILE_PATH
            output_file_path = os.path.join(tmpdir, 'output_file.mae')
            save_dir = tmpdir

            # Create the output file
            with open(output_file_path, 'w') as f:
                f.write('Test')

            with self.assertRaises(FileExistsError):
                convert_format(input_file_path, output_file_path, save_dir=save_dir)


    @patch('pyCADD.Dock.schrodinger.utils.StructureReader.read')
    def test_convert_format_unsupported_format(self, mock_read):
        # Test convert_format with an unsupported format
        with TemporaryDirectory() as tmpdir:
            input_file_path = TEST_PDB_FILE_PATH
            output_file_path = os.path.join(tmpdir, 'output_file.unkown')
            save_dir = tmpdir

            mock_read.side_effect = Exception

            with self.assertRaises(ValueError):
                convert_format(input_file_path, output_file_path, save_dir=save_dir)
            