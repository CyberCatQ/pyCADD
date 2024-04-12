import unittest

from pyCADD.Dock.schrodinger.config import (BaseConfig, DataConfig,
                                            DefaultDataConfig)


class TestBaseConfig(unittest.TestCase):
    def test_init(self):
        # Test initialization of BaseConfig
        base_config = BaseConfig()
        self.assertIsInstance(base_config, BaseConfig)


class TestDefaultDataConfig(unittest.TestCase):
    def test_init(self):
        # Test initialization of DefaultDataConfig with precision='SP'
        default_data_config = DefaultDataConfig()
        self.assertIsInstance(default_data_config, DefaultDataConfig)
        self.assertEqual(default_data_config.properties,
                         default_data_config.default_prop_SP)

        # Test initialization of DefaultDataConfig with precision='XP'
        default_data_config = DefaultDataConfig(precision='XP')
        self.assertIsInstance(default_data_config, DefaultDataConfig)
        self.assertEqual(default_data_config.properties,
                         default_data_config.default_prop_XP)

        # Test initialization of DefaultDataConfig with invalid precision
        with self.assertRaises(ValueError):
            DefaultDataConfig(precision='invalid')


class TestDataConfig(unittest.TestCase):
    def test_init(self):
        # Test initialization of DataConfig with precision='SP' and properties=None
        data_config = DataConfig(precision='SP')
        self.assertIsInstance(data_config, DataConfig)
        self.assertEqual(data_config.properties, data_config.default_prop_SP)

        # Test initialization of DataConfig with precision='XP' and properties=None
        data_config = DataConfig(precision='XP')
        self.assertIsInstance(data_config, DataConfig)
        self.assertEqual(data_config.properties, data_config.default_prop_XP)

        # Test initialization of DataConfig with precision='SP' and properties=list
        properties = ['property1', 'property2']
        data_config = DataConfig(precision='SP', properties=properties)
        self.assertIsInstance(data_config, DataConfig)
        self.assertEqual(data_config.properties,
                         data_config.default_prop_SP + properties)

        # Test initialization of DataConfig with precision='XP' and properties=list
        properties = ['property1', 'property2']
        data_config = DataConfig(precision='XP', properties=properties)
        self.assertIsInstance(data_config, DataConfig)
        self.assertEqual(data_config.properties,
                         data_config.default_prop_XP + properties)


if __name__ == '__main__':
    unittest.main()
