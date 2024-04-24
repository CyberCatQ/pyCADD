import os
import unittest
from tempfile import TemporaryDirectory

import pandas as pd

from pyCADD.Dock.schrodinger.common import MaestroFile
from pyCADD.Dock.schrodinger.ensemble import (_write_struc, get_docking_pairs,
                                              multi_dock, multi_extract_data,
                                              multi_grid_generate,
                                              multi_minimize, split_structure)
from pyCADD.utils.common import ChDir

from . import DEBUG, TEST_ASSETS_DIR


def get_test_pdbfiles(mappings, save_dir):
    pdbfiles = []
    with ChDir(save_dir):
        for mapping in mappings:
            pdb = mapping['pdb']
            ligand_name = mapping['ligand']
            pdbfile = MaestroFile(os.path.join(TEST_ASSETS_DIR, f"{pdb}.pdb"))
            ligand_searched = pdbfile.find_ligands([ligand_name])
            chain_id = ligand_searched[0].chain
            single_chain_st = pdbfile.get_chain_structure(chain_id=chain_id)
            single_chain_file = f"{pdb}-chain{chain_id}_{ligand_name}.pdb"
            single_chain_st.write(single_chain_file)
            pdbfiles.append(MaestroFile(single_chain_file))
    return pdbfiles


def get_mol_num_list(maestrofiles):
    mol_num_list = []
    for mae in maestrofiles:
        ligand = mae.find_ligands([mae.metadata.ligand_name])[0]
        mol_num_list.append(ligand.mol_num)
    return mol_num_list


def equal_file_list(file_list1, file_list2):
    file_list1 = sorted([file.file_path for file in file_list1])
    file_list2 = sorted(file_list2)
    return file_list1 == file_list2


class TestEnsemble(unittest.TestCase):
    def setUp(self) -> None:
        self.multi_ligands_file = os.path.join(
            TEST_ASSETS_DIR, 'structures', 'multiple_ligands.sdf')
        self.test_mappings = [
            {'receptor': 'P10275', 'pdb': '1XJ7', 'ligand': 'DHT'},
            {'receptor': 'P10275', 'pdb': '1XQ3', 'ligand': 'R18'},
            {'receptor': 'P10275', 'pdb': '2AM9', 'ligand': 'TES'},
            {'receptor': 'P10275', 'pdb': '2AM9', 'ligand': 'DTT'},
            {'receptor': 'P10275', 'pdb': '2YLP', 'ligand': 'TES'},
            {'receptor': 'P10275', 'pdb': '2YLP', 'ligand': '056'}
        ]

    def test_write_struc(self):
        with TemporaryDirectory() as tmpdir:
            save_dir = tmpdir if not DEBUG else os.path.join(
                os.getcwd(), 'tests')
            structure_index = 1
            structure_file = MaestroFile(self.multi_ligands_file)
            structure = structure_file.structures[0]
            result = _write_struc(structure_index,
                                  structure,
                                  file_prefix=structure_file.file_prefix,
                                  overwrite=True, save_dir=save_dir
                                  )
            expected_output_file = os.path.join(
                save_dir, f"{structure_file.file_prefix}-{structure_index}.maegz")
            self.assertEqual(result.file_path, expected_output_file)
            self.assertTrue(os.path.exists(expected_output_file))
            self.assertEqual(
                structure.property['i_user_StructureIndex'], structure_index)
            self.assertEqual(
                structure.property['s_user_StructureName'], f'{structure_file.file_prefix}-1')

    def test_split_structure(self):
        with TemporaryDirectory() as tmpdir:
            save_dir = tmpdir if not DEBUG else os.path.join(
                os.getcwd(), 'tests', 'multi_split')
            cpu_num = 4
            multi_structure_file = MaestroFile(self.multi_ligands_file)
            result = split_structure(
                multi_structure_file, save_dir=save_dir, overwrite=True, cpu_num=cpu_num)

            expected_output_files = [os.path.join(
                save_dir, f"{multi_structure_file.file_prefix}-{i}.maegz") for i in range(1, 4)]
            self.assertTrue(equal_file_list(result, expected_output_files))
            self.assertTrue(all([os.path.exists(f)
                            for f in expected_output_files]))

    def test_get_docking_pairs(self):
        grid_files = ['/path/to/grid1.zip', '/path/to/grid2.zip']
        ligand_files = ['/path/to/ligand1.mol2', '/path/to/ligand2.mol2']

        result = get_docking_pairs(grid_files, ligand_files)

        expected_result = [
            ('/path/to/grid1.zip', '/path/to/ligand1.mol2'),
            ('/path/to/grid1.zip', '/path/to/ligand2.mol2'),
            ('/path/to/grid2.zip', '/path/to/ligand1.mol2'),
            ('/path/to/grid2.zip', '/path/to/ligand2.mol2')
        ]
        self.assertEqual(result, expected_result)

    def test_multi_dock_workflow(self):
        with TemporaryDirectory() as tmpdir:

            # minimize
            ph = 7.4
            force_field = 'OPLS4'
            fill_side_chain = True
            add_missing_loop = True
            del_water = True
            watdist = 5.0
            rmsd_cutoff = 0.3
            overwrite = False
            cpu_num = 6
            save_dir = tmpdir if not DEBUG else os.path.join(
                os.getcwd(), 'tests', 'multi_minimize')
            structure_files = get_test_pdbfiles(self.test_mappings, save_dir)
            minimized_result = multi_minimize(
                structure_files,
                ph=ph,
                force_field=force_field,
                fill_side_chain=fill_side_chain,
                add_missing_loop=add_missing_loop,
                del_water=del_water,
                watdist=watdist,
                rmsd_cutoff=rmsd_cutoff,
                save_dir=save_dir,
                overwrite=overwrite,
                cpu_num=cpu_num
            )
            miminize_expected_output_files = [
                os.path.join(save_dir, f"{st_file.file_prefix}_minimized.mae") for st_file in structure_files]
            self.assertTrue(equal_file_list(minimized_result,
                            miminize_expected_output_files))
            self.assertTrue(all([os.path.exists(f)
                            for f in miminize_expected_output_files]))

            # generate grid
            box_center_molnum_list = get_mol_num_list(minimized_result)
            box_size_list = [20 for _ in range(6)]
            force_field = 'OPLS4'
            save_dir = tmpdir if not DEBUG else os.path.join(
                os.getcwd(), 'tests', 'multi_grid')
            overwrite = False
            cpu_num = 6
            grid_result = multi_grid_generate(
                minimized_result,
                box_center_molnum_list=box_center_molnum_list,
                box_size_list=box_size_list,
                force_field=force_field,
                save_dir=save_dir,
                overwrite=overwrite,
                cpu_num=cpu_num
            )
            grid_expected_output_files = [
                os.path.join(save_dir, f"{st_file.file_prefix}_glide-grid.zip") for st_file in structure_files]
            self.assertTrue(equal_file_list(
                grid_result, grid_expected_output_files))
            self.assertTrue(all([os.path.exists(f)
                            for f in grid_expected_output_files]))

            # dock
            save_dir = tmpdir if not DEBUG else os.path.join(
                os.getcwd(), 'tests', 'multi_dock')
            force_field = 'OPLS4'
            precision = 'SP'
            calc_rmsd = True
            include_receptor = False
            overwrite = False
            cpu_num = 6
            ligands = split_structure(
                self.multi_ligands_file, save_dir=save_dir, overwrite=overwrite, cpu_num=cpu_num)
            docking_pairs = [(grid, ligand)
                             for grid in grid_result for ligand in ligands]
            dock_result = multi_dock(
                docking_pairs,
                force_field=force_field,
                precision=precision,
                calc_rmsd=calc_rmsd,
                include_receptor=include_receptor,
                save_dir=save_dir,
                overwrite=overwrite,
                cpu_num=cpu_num
            )
            self.assertEqual(len(dock_result), 9)

            # data extract
            dock_data = multi_extract_data(dock_result, cpu_num=3)
            df = pd.DataFrame(dock_data)
            self.assertEqual(df.shape, (9, 19))


if __name__ == '__main__':
    unittest.main()
