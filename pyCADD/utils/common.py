import os

class BaseFile:
    '''
    基本文件类型

    Attributes
    ----------
    file_path : str
        文件绝对路径
    file_name : str
        纯文件名
    file_dir : str
        文件所在目录
    file_ext : str
        文件扩展名(不含有点)
    file_prefix : str
        文件前缀名
    file_suffix : str
        文件后缀名(与扩展名相同)
    '''
    def __init__(self, path) -> None:
        '''
        初始化
        
        Parameters
        ----------
        path : str
            文件路径
        '''

        if not os.path.exists(path):
            raise FileNotFoundError('File %s not found' % path)

        self.file_path = os.path.abspath(path)
        self.file_name = os.path.split(self.file_path)[-1]
        self.file_dir = os.path.split(self.file_path)[0]
        self.file_ext = os.path.splitext(self.file_name)[-1]
        self.file_prefix = os.path.splitext(self.file_name)[0]
        self.file_suffix = self.file_ext