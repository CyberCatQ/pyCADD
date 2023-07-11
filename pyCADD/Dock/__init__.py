import os
try:
    import schrodinger
except ImportError:
    print("Current environment does not include Schrodinger package. \nPlease use source to activate Schrodinger python environment first.")
    os._exit(1)

from pyCADD.Dock.console import Docker
from pyCADD.Dock.ensemble import _Console as MultiDocker