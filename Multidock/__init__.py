try:
    import pyCADD
except ImportError:
    import os
    os.system('python3 -m pip install pyCADD')
    
from pyCADD.utils.tool import init_log


init_log('pyCADD.Multidock')