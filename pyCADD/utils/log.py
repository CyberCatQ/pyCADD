import logging
import os
from datetime import datetime

from concurrent_log_handler import ConcurrentRotatingFileHandler
from rich.logging import RichHandler

def get_logfile_name():
    '''
    依据当前日期生成log文件的文件名
    如已有重复文件则添加递增后缀
    Return
    ---------
    str
        log文件名
    '''
    # 默认储存log文件的目录PATH
    log_dir = os.getcwd() + '/logs/'
    os.makedirs(log_dir, exist_ok=True)
    date = datetime.now()
    year = str(date.year)
    month = str(date.month)
    day = str(date.day)
    now = year + month.rjust(2, '0') + day.rjust(2, '0')

    i = 1
    while True:
        logfile = log_dir + now + '_' + str(i) + '.log'
        if os.path.exists(logfile):
            i += 1
        else:
            return logfile

def _init_log(logname):
    '''
    初始化logging配置
    '''
    logger = logging.getLogger(logname)
    logger.setLevel(logging.DEBUG)
    file_fmt = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    console_fmt = logging.Formatter('%(message)s')
    logfile = get_logfile_name()

    file_handler = ConcurrentRotatingFileHandler(logfile, 'a')
    file_handler.setFormatter(file_fmt)
    file_handler.setLevel(logging.DEBUG)

    console_handler = RichHandler()
    console_handler.setFormatter(console_fmt)
    console_handler.setLevel(logging.INFO)

    logger.addHandler(file_handler)
    logger.addHandler(console_handler)
    logger.logfilename = logfile

    return logger
