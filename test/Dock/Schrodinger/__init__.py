import logging
import os

TEST_ASSETS_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), 'assets')
TEST_PDB_FILE_PATH = os.path.join(TEST_ASSETS_DIR, '3OAP.pdb')
DEBUG = os.environ.get('PYCADD_DEBUG', False)

def init_logger():
    # reset logger during test
    logger = logging.getLogger('pyCADD')
    for handler in logger.handlers:
        logger.removeHandler(handler)
    handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter('%(message)s'))
    logger.setLevel(logging.WARNING)
    logger.addHandler(handler)
