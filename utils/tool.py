import os 
import logging
logger = logging.getLogger('pyCADD.utils.tool')

def mkdirs(path_list:list):
    '''
    获取包含多个PATH的列表 尝试创建列表中的所有目录
    
    Parameter
    ----------
    path_list : list
        要创建的PATH列表
    '''

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
    log_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/logs/'

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