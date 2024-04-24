import os
import unittest

from pyCADD.Dock.common import (EnsembleInputFile, PDBFile, PDBLine,
                                PDBLineParser)

from . import TEST_ASSETS_DIR, TEST_PDB_FILE_PATH

TEST_PDB_LINE_NUM = 2373


class TestPDBLine(unittest.TestCase):
    def test_init(self):
        line = "ATOM      1  N   ASP A   1      10.000  20.000  30.000  1.00  0.00           N"
        pdb_line = PDBLine(line)
        self.assertEqual(pdb_line.record_name, "ATOM")
        self.assertEqual(pdb_line.atom_idx, "1")
        self.assertEqual(pdb_line.atom_name, "N")
        self.assertEqual(pdb_line.alt_loc, "")
        self.assertEqual(pdb_line.res_name, "ASP")
        self.assertEqual(pdb_line.chain_id, "A")
        self.assertEqual(pdb_line.res_id, "1")
        self.assertEqual(pdb_line.insertion_code, "")
        self.assertEqual(pdb_line.coord_x, "10.000")
        self.assertEqual(pdb_line.coord_y, "20.000")
        self.assertEqual(pdb_line.coord_z, "30.000")
        self.assertEqual(pdb_line.occupancy, "1.00")
        self.assertEqual(pdb_line.temp_factor, "0.00")
        self.assertEqual(pdb_line.element, "N")
        self.assertEqual(pdb_line.charge, "")

    def test_get_atom_name(self):
        line = "ATOM      1  1HB ASP A   1      10.000  20.000  30.000  1.00  0.00           N"
        pdb_line = PDBLine(line)
        self.assertEqual(pdb_line.get_atom_name(), "HB1")


class TestPDBLineParser(unittest.TestCase):
    def test_init_with_pdb_str(self):
        pdb_str = "ATOM      1  N   ASP A   1      10.000  20.000  30.000  1.00  0.00           N\nATOM      2  CA  ASP A   1      11.000  21.000  31.000  1.00  0.00           C"
        pdb_parser = PDBLineParser(pdb_str=pdb_str)
        self.assertEqual(len(pdb_parser.get_lines()), 2)

    def test_init_with_pdb_file(self):
        pdb_file = TEST_PDB_FILE_PATH
        pdb_parser = PDBLineParser(pdb_file=pdb_file)
        self.assertEqual(len(pdb_parser.get_lines()), TEST_PDB_LINE_NUM)

class TestPDBFile(unittest.TestCase):
    def test_init(self):
        pdb_file = PDBFile(TEST_PDB_FILE_PATH)
        self.assertEqual(pdb_file.pdbid, "3OAP")
        self.assertEqual(len(pdb_file.pdb_parser.get_lines()), TEST_PDB_LINE_NUM)

    def test_get_lines(self):
        pdb_file = PDBFile(TEST_PDB_FILE_PATH)
        lines = pdb_file.get_lines()
        self.assertIsInstance(lines, list)
        self.assertEqual(len(lines), TEST_PDB_LINE_NUM)

    def test_get_chain(self):
        pdb_file = PDBFile(TEST_PDB_FILE_PATH)
        chain = pdb_file.get_chain("A", return_str=True)
        self.assertIsInstance(chain, str)
        chain = pdb_file.get_chain("A")
        self.assertIsInstance(chain, list)
        self.assertEqual(chain[0], "ATOM      1  N   GLU A 228      49.559  71.385  12.045  1.00 64.80           N  ")
        
class TestEnsembleInputFile(unittest.TestCase):
    def setUp(self):
        self.csv_file_path = os.path.join(TEST_ASSETS_DIR, 'ensemble', 'P10275.csv')
        self.ini_file_path = os.path.join(TEST_ASSETS_DIR, 'ensemble', 'P10275.ini')
        self.yaml_file_path = os.path.join(TEST_ASSETS_DIR, 'ensemble','P10275.yml')
        self.expected_mapping = [{'receptor': 'P10275', 'pdb': '1XJ7', 'ligand': 'DHT'},
                                  {'receptor': 'P10275', 'pdb': '1XQ3', 'ligand': 'R18'},
                                  {'receptor': 'P10275', 'pdb': '2AM9', 'ligand': 'TES'},
                                  {'receptor': 'P10275', 'pdb': '2AM9', 'ligand': 'DTT'},
                                  {'receptor': 'P10275', 'pdb': '2YLP', 'ligand': 'TES'},
                                  {'receptor': 'P10275', 'pdb': '2YLP', 'ligand': '056'}]

    def test_from_csv(self):
        ensemble_file = EnsembleInputFile.from_csv(self.csv_file_path, header=True)
        self.assertEqual(ensemble_file.file_path, self.csv_file_path)
        self.assertEqual(ensemble_file.mappings, self.expected_mapping)

    def test_from_ini(self):
        ensemble_file = EnsembleInputFile.from_ini(self.ini_file_path)
        self.assertEqual(ensemble_file.file_path, self.ini_file_path)
        self.assertEqual(ensemble_file.mappings, self.expected_mapping)

    def test_from_yaml(self):
        ensemble_file = EnsembleInputFile.from_yaml(self.yaml_file_path)
        self.assertEqual(ensemble_file.file_path, self.yaml_file_path)
        self.assertEqual(ensemble_file.mappings, self.expected_mapping)

    def test_parse_file_csv(self):
        ensemble_file = EnsembleInputFile.parse_file(self.csv_file_path, header=True)
        self.assertEqual(ensemble_file.file_path, self.csv_file_path)
        self.assertEqual(ensemble_file.mappings, self.expected_mapping)

    def test_parse_file_ini(self):
        ensemble_file = EnsembleInputFile.parse_file(self.ini_file_path)
        self.assertEqual(ensemble_file.file_path, self.ini_file_path)
        self.assertEqual(ensemble_file.mappings, self.expected_mapping)

    def test_parse_file_yaml(self):
        ensemble_file = EnsembleInputFile.parse_file(self.yaml_file_path)
        self.assertEqual(ensemble_file.file_path, self.yaml_file_path)
        self.assertEqual(ensemble_file.mappings, self.expected_mapping)

    def test_get_pairs_list(self):
        ensemble_file = EnsembleInputFile(self.csv_file_path)
        ensemble_file.mappings = self.expected_mapping
        pairs_list = ensemble_file.get_pairs_list()
        expected_pairs_list = [('1XJ7', 'DHT'), ('1XQ3', 'R18'), ('2AM9', 'TES'),
                               ('2AM9', 'DTT'), ('2YLP', 'TES'), ('2YLP', '056')]
        self.assertTrue(all([pair in expected_pairs_list for pair in pairs_list]))
        self.assertEqual(len(pairs_list), 6)

    def test_get_pdbid_list(self):
        ensemble_file = EnsembleInputFile(self.csv_file_path)
        ensemble_file.mappings = self.expected_mapping
        pdbid_list = ensemble_file.get_pdbid_list()
        expected_pdbid_list = ['1XJ7', '1XQ3', '2AM9', '2YLP']
        self.assertEqual(pdbid_list, expected_pdbid_list)
        self.assertEqual(len(pdbid_list), 4)

    def test_get_ligand_list(self):
        ensemble_file = EnsembleInputFile(self.csv_file_path)
        ensemble_file.mappings = self.expected_mapping
        ligand_list = ensemble_file.get_ligand_list()
        expected_ligand_list = ['056', 'DHT', 'DTT', 'R18', 'TES']
        self.assertEqual(ligand_list, expected_ligand_list)
        self.assertEqual(len(ligand_list), 5)

if __name__ == '__main__':
    unittest.main()
