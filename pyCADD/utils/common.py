import os
from configparser import ConfigParser
from threading import Thread


class BaseFile:

    def __init__(self, path: str) -> None:
        """
        Basic file class

        Attributes:
            file_path (str): Absolute file path
            file_name (str): File basename without directory path
            file_dir (str): File directory path
            file_ext (str): File extension name without dot
            file_prefix (str): File prefix name
            file_suffix (str): File suffix name without dot

        Args:
            path (str): file path string
        """
        if isinstance(path, BaseFile):
            path = path.file_path
        self.file_path = os.path.abspath(path)
        self.file_name = os.path.basename(self.file_path)
        self.file_dir = os.path.split(self.file_path)[0]
        self.file_ext = os.path.splitext(self.file_name)[-1].replace(".", "")
        self.file_prefix = os.path.splitext(self.file_name)[0]
        self.file_suffix = self.file_ext


class File(BaseFile):
    def __init__(self, path: str, exist: bool = True) -> None:
        """File class

        Attributes:
            file_path (str): Absolute file path
            file_name (str): File basename without directory path
            file_dir (str): File directory path
            file_ext (str): File extension name without dot
            file_prefix (str): File prefix name
            file_suffix (str): File suffix name without dot

        Args:
            path (str): file path string
            exist (bool, optional): Check if the file exists. Defaults to True.

        Example:
            >>> f = File("dir_name/example.txt")
            >>> f.file_path
            '/home/user/dir_name/example.txt'
            >>> f.file_name
            'example.txt'
            >>> f.file_dir
            '/home/user/dir_name'
            >>> f.file_ext
            'txt'
            >>> f.file_prefix
            'example'
            >>> f.file_suffix
            'txt'

        Raises:
            FileNotFoundError: If the file does not exist and exist is True
        """
        super().__init__(path)
        if exist and not os.path.exists(self.file_path):
            raise FileNotFoundError(f"File not found: {self.file_path}")

    def __str__(self) -> str:
        return f"<File at {self.file_path} exist={os.path.exists(self.file_path)}>"

    def __repr__(self) -> str:
        return self.__str__()

    def read(self, binary: bool = False) -> str:
        """Read the file content

        Returns:
            str: File content string
        """
        with open(self.file_path, "rb" if binary else "r") as f:
            return f.read()


class FixedConfig(ConfigParser):
    def __init__(self, defaults=None):
        """Fixed optionxform method for ConfigParser due to upper string conversion"""
        ConfigParser.__init__(self, defaults=defaults)

    def optionxform(self, optionstr):
        return optionstr


class FixedThread(Thread):
    def __init__(self, *args, **kwargs):
        """Thread class with exception handling"""
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


class ChDir:
    def __init__(self, path: str = None, exist: bool = False, delete: bool = False):
        """Change working directory to the given path

        Args:
            path (str): path to change, None to create a temporary directory.
            exist (bool, optional): Whether the directory must exist. False will create the directory if it does not exist. Defaults to False.
            delete (bool, optional): Whether to delete the directory after use. Defaults to False.
        """
        from tempfile import mkdtemp

        self.path = path if path is not None else mkdtemp()
        self.delete = delete
        self.cwd = os.getcwd()
        if not os.path.exists(self.path):
            if exist:
                raise FileNotFoundError(f"Directory {self.path} not found")
            else:
                os.makedirs(self.path, exist_ok=True)

    def __enter__(self):
        os.chdir(self.path)

    def __exit__(self, exc_type, exc_value, traceback):
        if self.delete:
            import shutil

            shutil.rmtree(self.path)
        os.chdir(self.cwd)
