import os
import time
import logging
import requests
import multiprocessing

# For Schrodinger 2021-2 or newer release
import importlib
importlib.reload(multiprocessing)
from typing import Iterable

from rich.progress import SpinnerColumn, TextColumn, BarColumn, Progress, TimeElapsedColumn, TimeRemainingColumn
from rich.table import Column
from configparser import ConfigParser
from threading import Thread
from concurrent.futures import ProcessPoolExecutor, Future

NUM_PARALLEL = multiprocessing.cpu_count() // 4 * 3
logger = logging.getLogger(__name__)

def _multiprocssing_run(func, _iterable:Iterable, *args, job_name:str, num_parallel:int):
    '''
    多进程运行函数 并自动生成进度条

    Parameters
    ----------
    func : function
        运行函数
    _iterable : Iterable
        可迭代对象 即函数作用对象总集
    *args
        传入func函数的其他参数
    job_name : str
        进程名称
    num_parallel : int
        进程数量
    
    Returns
    -------
    List
        所有成功完成任务的返回值
    
    '''
    progress, taskID = _get_progress(job_name, 'bold cyan', len(_iterable))
    returns = []
    def success_handler(future:Future):
        '''
        处理完成的任务(成功/失败)
        '''
        if future.exception() is not None:
            logger.debug(f'Multiprocessing Warnning: {future.exception()}')
        else:
            returns.append(future.result())
        progress.update(taskID, advance=1)
    '''
    def _error_handler(error:Exception):
        
        异常处理函数
        
        logger.debug(f'Warnning: {error}')
        progress.update(taskID, advance=1)
    '''
    progress.start()
    progress.start_task(taskID)

    with ProcessPoolExecutor(max_workers=num_parallel) as pool:
        for item in _iterable:
            if isinstance(item, Iterable):
                future = pool.submit(func, *item, *args)
                future.add_done_callback(success_handler)
            else:
                future = pool.submit(func, item, *args)
                future.add_done_callback(success_handler)

    # 当ligands数量较多时 父进程将卡死于pool.join()
    # bug原因未知
    # pool = multiprocessing.Pool(num_parallel, maxtasksperchild=1)
    # for item in _iterable:
    #    if isinstance(item, Iterable):
    #        pool.apply_async(func, args=(*item, *args), callback=success_handler, error_callback=_error_handler)
    #    else:
    #        pool.apply_async(func, args=(item, *args), callback=success_handler, error_callback=_error_handler)
    # pool.close()
    # pool.join()

    progress.stop()
    return returns
 
def get_lib_dir():
    '''
    获取pyCADD库所在的Absolute PATH
    '''
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def download_pdb(pdbid, download_dir:str=None, overwrite:bool=False) -> None:
    '''
    从RCSB服务器下载PDB文件

    Parameters
    ----------
    pdbid : str
        PDB ID
    download_dir : str
        下载目录
    overwrite : bool
        是否覆盖已存在的文件

    '''
    base_url = 'https://files.rcsb.org/download/'
    pdbfile = pdbid + '.pdb'
    download_dir = os.getcwd() if download_dir is None else download_dir
    downloaded_file = os.path.join(download_dir, pdbfile)

    if os.path.exists(downloaded_file) and not overwrite:
        return

    logger.debug('Downloading %s ...' % pdbid)
    url = base_url + pdbid + '.pdb'
    response = requests.get(url)
    pdb_data = response.text
    with open(downloaded_file, 'w') as f:
        f.write(pdb_data)
    
    logger.debug('%s.pdb downloaded.' % pdbid)

def download_pdb_list(pdblist:list, download_dir:str=None, overwrite:bool=False) -> None:
    '''
    多线程下载PDB ID列表中的所有PDB文件
    Parameters
    ----------
    pdblist : list
        PDB列表
    download_dir : str
        下载目录
    overwrite : bool
        是否覆盖已存在的文件
    '''
    download_dir = os.getcwd() if download_dir is None else download_dir
    threads = []
    for pdbid in pdblist:
        t = Thread(target=download_pdb, args=(pdbid, download_dir, overwrite))
        threads.append(t)
    for t in threads:
        t.start()
    for t in threads:
        t.join()

def makedirs_from_list(path_list: list) -> None:
    '''
    输入包含多个PATH的列表 尝试创建列表中的所有目录

    Parameters
    ----------
    path_list : list
        要创建的PATH列表
    '''
    for path in path_list:
        os.makedirs(path, exist_ok=True)

def _get_progress(name: str, description: str, total: int, start:bool=False):
    '''
    创建Rich进度条

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

# 废弃
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

def get_config(config_file: str) -> Myconfig:
    '''
    读取配置文件

    Parameters
    ----------
    config_file : str
        配置文件路径

    Return
    ----------
    Myconfig
        配置文件对象
    '''
    config = Myconfig()
    config.read(config_file)
    return config