'''
Data Analyzer for Computer-aid drug design
'''
from pandas import DataFrame
from pyCADD.Dance.core import _Evaluator

class Evaluator(_Evaluator):
    '''
    Evaluator class for pyCADD.Dance.
    '''
    def __init__(self, train_data:DataFrame=None, test_data:DataFrame=None, data:DataFrame=None, test_size=0.25, random_state=42) -> None:
        if train_data is not None and test_data is not None:
            super().__init__(train_data, test_data)
        else:
            if data is not None:
                train_data, test_data = self._preprocess_data(data, test_size, random_state)
                super().__init__(train_data, test_data)
            else:
                raise ValueError('Either train_data/test_data or data must be provided.')