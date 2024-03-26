import os
from pyCADD.utils.tool import is_gaussian_available
if not is_gaussian_available():
    print('Gaussian 16 is not installed or not in PATH.\nDensity module may not work properly.')