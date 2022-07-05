#!/bin/sh
# Conversion of .mae files to .pdb files
$SCHRODINGER/run python3 -c "from pyCADD.Dock.common import MaestroFile; MaestroFile.convert_format('$1', 'pdb')"