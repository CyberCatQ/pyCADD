
import os
base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
log_dir = base_dir + '/logs'

# 配置log
import logging
logger = logging.getLogger('pyCADD')
logger.setLevel(level = logging.INFO)
file_fmt = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
console_fmt = logging.Formatter('%(levelname)s - %(message)s')

filehandler = logging.FileHandler(log_dir + '/%s.log' % now, 'a')
filehandler.setLevel(logging.INFO)
filehandler.setFormatter(file_fmt)
consolehandler = logging.StreamHandler()
consolehandler.setLevel(logging.INFO)
consolehandler.setFormatter(console_fmt)

logger.addHandler(filehandler)
logger.addHandler(consolehandler)

from pyCADD.Multidock import core

class Multidock:
    '''
    Multi Mode
    '''

    def __init__(self) -> None:
        self.list_file = ''         # 单一基因的全部PDB晶体列表文件
        self.ligand_file = ''       # 外源配体文件
        self.precision = ''         # 对接精度
        self.cpus = ''              # 对接工作使用的进程数
        self.flag = ''              # 对接工作启动确认标记
        self.mmgbsaFlag = ''        # 是否计算MM-GBSA结合能的标记
        self.list_filename = ''     # 基因名称(晶体列表文件的名称)
        self.notpass = []           # 内源配体对接不成功的晶体列表
        self.dock_fail = []         # 外源配体对接不成功的晶体列表

    
    