import os
import shutil
import unittest
from tempfile import TemporaryDirectory
from unittest.mock import MagicMock, Mock, patch

from pyCADD.Dock.schrodinger.common import MaestroFile
from pyCADD.Dock.schrodinger.core import dock, grid_generate, minimize
from pyCADD.Dock.schrodinger.data import (extract_docking_data,
                                          save_docking_data)
from pyCADD.utils.common import ChDir

from . import DEBUG, TEST_PDB_FILE_PATH, init_logger

init_logger()


class TestData(unittest.TestCase):

    def _dock(sel, include_recep: bool = False):
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

    def test_extract_docking_data_with_recep(self):
        # Test extracting docking data from structure with receptor
        with TemporaryDirectory() as tmpdir:
            testdir = os.path.join(os.getcwd(), 'tests') if DEBUG else tmpdir
            with ChDir(testdir):
                docking_result = self._dock(True)
            result_data_list = extract_docking_data(docking_result)
            self.assertEqual(len(result_data_list), 1)
            self.assertEqual(
                result_data_list[0]['pdbid'], docking_result.metadata.pdbid)
            self.assertEqual(
                result_data_list[0]['precision'], docking_result.metadata.precision)
            self.assertEqual(
                result_data_list[0]['internal_ligand_name'], docking_result.metadata.internal_ligand_name)
            self.assertEqual(
                result_data_list[0]['docking_ligand_name'], docking_result.metadata.docking_ligand_name)
            self.assertLess(result_data_list[0]['Docking_Score'], -8.0)

    def test_save_docking_data(self):
        # Test saving docking data to a CSV file from the structure without receptor
        with TemporaryDirectory() as tmpdir:
            testdir = os.path.join(os.getcwd(), 'tests') if DEBUG else tmpdir

            with ChDir(testdir):
                docking_result = self._dock()

            data_list = docking_result.get_results()
            df_mock = self.create_mock_dataframe(data_list)

            with patch('pandas.DataFrame', return_value=df_mock) as mock_dataframe:
                save_dir = testdir
                save_docking_data(
                    docking_result, save_dir=save_dir, overwrite=True)
                output_file = os.path.join(
                    save_dir, f"{docking_result.file_prefix}.csv")
                self.assertTrue(os.path.exists(output_file))
                mock_dataframe.assert_called_once_with(data_list)
                df_mock.to_csv.assert_called_once_with(
                    output_file, index=False)
                with self.assertRaises(FileExistsError):
                    save_docking_data(docking_result, save_dir=save_dir)

    def create_mock_dataframe(self, data_list):
        mock_dataframe = MagicMock()
        mock_dataframe.return_value = mock_dataframe
        mock_dataframe.to_csv = Mock()
        mock_dataframe.to_csv.side_effect = mock_to_csv
        mock_dataframe.__getitem__.side_effect = lambda key: data_list[key]
        return mock_dataframe
    
def mock_to_csv(output_file, index=False):
    with open(output_file, 'w') as f:
        f.write('mock data')

if __name__ == '__main__':
    unittest.main()
