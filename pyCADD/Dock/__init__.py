import sys
try:
    import schrodinger
except ImportError:
    print("Current environment does not include Schrodinger package. \nPlease activate Schrodinger python environment first.")
    sys.exit(1)

from pyCADD.Dock.console import Docker
from pyCADD.Dock.ensemble import _Console as MultiDocker