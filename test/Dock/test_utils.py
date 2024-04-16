import unittest
from unittest.mock import patch

from pyCADD.Dock.utils import check_pdb, get_input_pdbid


class TestUtils(unittest.TestCase):
    def test_check_pdb_valid(self):
        # Test check_pdb function with a valid PDB ID
        pdb_id = '3OAP'
        result = check_pdb(pdb_id)
        self.assertTrue(result)

    def test_check_pdb_invalid(self):
        # Test check_pdb function with an invalid PDB ID
        pdb_id = 'invalid_pdb_id'
        result = check_pdb(pdb_id)
        self.assertFalse(result)

    def test_get_input_pdbid_from_directory(self):
        # Test get_input_pdbid function when PDB ID is obtained from directory name
        with patch('os.getcwd', return_value='/path/to/3OAP'):
            result = get_input_pdbid()
        self.assertEqual(result, '3OAP')

    @patch('builtins.input', return_value='3OAP')
    def test_get_input_pdbid_from_user_input(self, mock_input):
        # Test get_input_pdbid function when PDB ID is obtained from user input
        result = get_input_pdbid()
        self.assertEqual(result, '3OAP')