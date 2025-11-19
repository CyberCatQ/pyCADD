import logging
import os

from .. import TEST_ASSETS_DIR, TEST_PDB_FILE_PATH, TEST_HETATM_NAME, init_logger

DEBUG = os.environ.get('PYCADD_DEBUG', False)
