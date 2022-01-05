import os
import logging
import time

from cloghandler import ConcurrentRotatingFileHandler
from rich.logging import RichHandler
from rich.progress import SpinnerColumn, TextColumn, BarColumn, Progress, TimeElapsedColumn, TimeRemainingColumn
from rich.table import Column
from configparser import ConfigParser


def mkdirs(path_list: list):
    '''
    获取包含多个PATH的列表 尝试创建列表中的所有目录

    Parameter
    ----------
    path_list : list
        要创建的PATH列表
    '''
    logger = logging.getLogger('pyCADD.utils.tool')
    for path in path_list:
        try:
            os.mkdir(path)
            logger.info('Created directory %s' % path)
        except FileExistsError:
            continue


def generate_logfile_name():
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
    mkdirs([log_dir])

    from datetime import datetime
    # 获取今日日期
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
    初始化配置log
    '''
    logger = logging.getLogger(logname)
    logger.setLevel(level=logging.DEBUG)
    file_fmt = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    console_fmt = logging.Formatter('%(message)s')

    logfile = generate_logfile_name()
    filehandler = ConcurrentRotatingFileHandler(logfile, 'a')
    filehandler.setLevel(logging.DEBUG)
    filehandler.setFormatter(file_fmt)
    # consolehandler = logging.StreamHandler()
    consolehandler = RichHandler()
    consolehandler.setLevel(logging.INFO)
    consolehandler.setFormatter(console_fmt)

    logger.addHandler(filehandler)
    logger.addHandler(consolehandler)
    logger.logfilename = logfile

    return logger

def init_log(logname):
    '''
    getLogger
    '''
    return logging.getLogger(logname)


def _get_progress(name: str, description: str, total: int, start:bool=False):
    '''
    创建进度条

    Parameters
    ----------
    name : str
        进度条进程名
    description : str
        进度条名称样式(example: 'bold red')
    total : int
        标志项目总进度为100%时的长度

    Reture
    ----------
    rich.progress.Progress, task ID
        进度条对象, 任务ID(用于update)
    '''

    text_column = TextColumn("{task.description}",
                             table_column=Column(), justify='right')
    percent_column = TextColumn(
        "[bold green]{task.percentage:.1f}%", table_column=Column())
    finished_column = TextColumn(
        "[bold purple]{task.completed} of {task.total}")
    bar_column = BarColumn(bar_width=None, table_column=Column())
    progress = Progress(SpinnerColumn(), text_column, "•", TimeElapsedColumn(
    ), "•", percent_column, bar_column, finished_column, TimeRemainingColumn())

    task = progress.add_task('[%s]%s' % (
        description, name), total=total, start=start)

    return progress, task


def check_file_update_progress(file_path: str, progress:Progress, task_ID: str, time_sleep: int = 3):
    '''
    定时检查文件是否存在 已存在则更新进度条

    Parameters
    ----------
    file_path : str
        检查的文件路径
    progress : rich.progress.Progress
        进度条对象
    task_ID : str
        要更新的进度条任务ID
    time_sleep : int
        检查间隔时间

    '''
    while os.path.exists(file_path) == False:
        time.sleep(time_sleep)
    
    progress.update(task_ID, advance=1)
    time.sleep(0.5)


class Myconfig(ConfigParser):
    '''
    重写以解决配置读取大小写修改问题
    '''

    def __init__(self, defaults=None):
        ConfigParser.__init__(self, defaults=defaults)

    def optionxform(self, optionstr):
        return optionstr


def tail_progress(file:str, callback, wait=0.1):
    '''
    追踪文件末尾行内容

    Parameters
    ----------
    file : str
        追踪的文件
    callback : func
        回调函数 函数返回True时结束追踪
    wait : int | float
        追踪间隔时长(s)
    '''

    with open(file) as f:
        f.seek(0, os.SEEK_END)
        while True:
            line = f.readline()
            if line:
                if callback(line):
                    break
            time.sleep(wait)