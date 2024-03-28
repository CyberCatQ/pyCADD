import os
from configparser import ConfigParser
from threading import Thread

class BaseFile:
    """
    Basic file class

    Attributes:
        file_path (str): Absolute file path
        file_name (str): File basename without directory path
        file_dir (str): File directory path
        file_ext (str): File extension name without dot
        file_prefix (str): File prefix name
        file_suffix (str): File suffix name without dot
    """

    def __init__(self, path: str, exist: bool = True) -> None:
        """Initialize file object

        Args:
            path (str): file path string
            exist (bool, optional): Check if the file exists. Defaults to True.

        Raises:
            FileNotFoundError: If the file does not exist and exist is True
        """
        if exist and not os.path.exists(path):
            raise FileNotFoundError('File %s not found' % path)

        self.file_path = os.path.abspath(path)
        self.file_name = os.path.basename(self.file_path)
        self.file_dir = os.path.split(self.file_path)[0]
        self.file_ext = os.path.splitext(self.file_name)[-1].replace('.', '')
        self.file_prefix = os.path.splitext(self.file_name)[0]
        self.file_suffix = self.file_ext


class File(BaseException):
    def __init__(self, *args: object) -> None:
        super().__init__(*args)


class FixedConfig(ConfigParser):
    """Fixed optionxform method for ConfigParser due to upper string conversion"""

    def __init__(self, defaults=None):
        ConfigParser.__init__(self, defaults=defaults)

    def optionxform(self, optionstr):
        return optionstr
    
class FixedThread(Thread):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._exception = None
    
    def run(self):
        try:
            if self._target:
                self._target(*self._args, **self._kwargs)
        except Exception as e:
            self._exception = e
    
    def join(self):
        super().join()
        if self._exception:
            raise self._exception

class TimeoutError(Exception):
    pass