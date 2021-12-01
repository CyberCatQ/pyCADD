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
    