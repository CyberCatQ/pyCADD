import unittest

from pyCADD.Dock.schrodinger.config import (DataConfig, DefaultDataConfig,
                                            SPConfig, XPConfig)


class TestDefaultDataConfig(unittest.TestCase):
    def test_init(self):
        # Test initialization of DefaultDataConfig with precision='SP'
        default_data_config = DefaultDataConfig()
        self.assertIsInstance(default_data_config, DefaultDataConfig)


class TestDataConfig(unittest.TestCase):
    def test_init(self):
        # Test initialization of DataConfig with precision='SP' and properties=None
        data_config = DataConfig(precision='SP')
        self.assertIsInstance(data_config, DataConfig)
        self.assertEqual(data_config.properties, SPConfig().properties)

        # Test initialization of DataConfig with precision='XP' and properties=None
        data_config = DataConfig(precision='XP')
        self.assertIsInstance(data_config, DataConfig)
        self.assertEqual(data_config.properties, XPConfig().properties)

        # Test initialization of DataConfig with precision=None and properties=None
        data_config = DataConfig()
        self.assertIsInstance(data_config, DataConfig)
        self.assertEqual(data_config.properties,
                         DefaultDataConfig().properties)

        # Test initialization of DataConfig with precision=None and additional properties
        data_config = DataConfig(properties={'new_property': 'new_value'})
        self.assertIsInstance(data_config, DataConfig)
        test_dict = DefaultDataConfig().properties.copy()
        test_dict.update({'new_property': 'new_value'})
        self.assertEqual(data_config.properties, test_dict)

    def test_items(self):
        # Test items() method
        data_config = DataConfig(precision='SP')
        self.assertEqual(data_config.items(), SPConfig().properties.items())

    def test_keys(self):
        # Test keys() method
        data_config = DataConfig(precision='SP')
        self.assertEqual([k for k in data_config.keys()], [
                         k for k in SPConfig().properties.keys()])

    def test_values(self):
        # Test values() method
        data_config = DataConfig(precision='SP')
        self.assertEqual([v for v in data_config.values()], [
                         v for v in SPConfig().properties.values()])

    def test_get(self):
        # Test get() method
        data_config = DataConfig(precision='SP')
        self.assertEqual(data_config.get('Docking_Score'),
                         SPConfig().properties.get('Docking_Score'))

    def test_getitem(self):
        # Test __getitem__() method
        data_config = DataConfig(precision='SP')
        self.assertEqual(data_config['Docking_Score'],
                         SPConfig().properties['Docking_Score'])

    def test_setitem(self):
        # Test __setitem__() method
        data_config = DataConfig(precision='SP')
        data_config['Docking_Score'] = 'new_value'
        self.assertEqual(data_config['Docking_Score'], 'new_value')


class TestSPConfig(unittest.TestCase):
    def test_properties(self):
        sp_config = SPConfig()
        self.assertEqual(sp_config.properties, {
            'Docking_Score': 'r_i_docking_score',
            'rmsd': 'r_i_glide_rmsd_to_input',
            'ligand_efficiency': 'r_i_glide_ligand_efficiency',
            'rotatable_bonds': 'i_i_glide_rotatable_bonds',
            'ecoul': 'r_i_glide_ecoul',
            'evdw': 'r_i_glide_evdw',
            'emodel': 'r_i_glide_emodel',
            'energy': 'r_i_glide_energy',
            'einternal': 'r_i_glide_einternal',
            'lipo': 'r_i_glide_lipo',
            'hbond': 'r_i_glide_hbond',
            'metal': 'r_i_glide_metal',
            'rewards': 'r_i_glide_rewards',
            'erotb': 'r_i_glide_erotb',
            'esite': 'r_i_glide_esite'
        })


class TestXPConfig(unittest.TestCase):
    def test_properties(self):
        xp_config = XPConfig()
        self.assertEqual(xp_config.properties, {
            'Docking_Score': 'r_i_docking_score',
            'rmsd': 'r_i_glide_rmsd_to_input',
            'ligand_efficiency': 'r_i_glide_ligand_efficiency',
            'rotatable_bonds': 'i_i_glide_rotatable_bonds',
            'ecoul': 'r_i_glide_ecoul',
            'evdw': 'r_i_glide_evdw',
            'emodel': 'r_i_glide_emodel',
            'energy': 'r_i_glide_energy',
            'einternal': 'r_i_glide_einternal',
            'XP_Hbond': 'r_glide_XP_HBond'
        })


if __name__ == '__main__':
    unittest.main()
