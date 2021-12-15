import logging

from pyCADD.ui import UI
try:
    from pyCADD.VSW.base import VSW
except ImportError:
    import os
    os.system('run python3 -m pip install rich ConcurrentLogHandler')
    from pyCADD.VSW.base import VSW

logger = logging.getLogger('pyCADD.VSW')

class VSW_UI(UI):
    '''
    虚拟筛选UI
    '''
    def __init__(self, menu_name: str = 'VSW') -> None:
        super().__init__(menu_name=menu_name)
    

    