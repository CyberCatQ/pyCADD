import os
import platform
from datetime import datetime

from rich import box, print
from rich.console import Group
from rich.padding import Padding
from rich.panel import Panel
from rich.table import Column, Table
from rich.text import Text
from rich.prompt import Confirm, Prompt

date = datetime.now()
year = str(date.year)
month = str(date.month).rjust(2,'0')
day = str(date.day).rjust(2,'0')
now = "%s-%s-%s" % (year, month, day)

__version__= "Undefined"
__update_date__ = "Undefined"

for line in open(os.path.dirname(os.path.abspath(__file__)) + '/__init__.py'):
    if line.startswith('__version__') or line.startswith("__update_date__"):
        exec(line.strip())
    
class UI:
    '''
    pyCADD程序用户交互界面(user interface)
    '''

    def __init__(self, menu_name: str = 'Main') -> None:
        self.version = __version__
        self.update_date = __update_date__
        self.menu_name = '[bold magenta]Menu: %s' % menu_name
        self.options = ''
        self.additional_info = ''
        self.schrodinger = Text(os.environ['SCHRODINGER'], style='u i')
        self.gauss_dir = Text(os.path.dirname(os.popen('which g16').read()), style='u i')

        self._additional_info_dict = {}
        self._info_index = 0

        self.gauss_check = True
        self.schrodinger_check = True
        if not self.schrodinger:
            self.schrodinger = Text('Not installed', style='bold red')
            self.schrodinger_check = False
        
        if not self.gauss_dir:
            self.gauss_dir = Text('Not installed', style='bold red')
            self.gauss_check = False
            
    @property
    def title(self) -> None:
        '''
        程序标题样式
        '''
        return Text.assemble(
            ('pyCADD', 'bold medium_spring_green'), 
            ' -- A ', ('Python Package', 'bold yellow'), 
            ' For ', 
            ('Computer-aid Drug Design', 'bold cyan')
            )

    @property
    def basic_info(self) -> None:
        '''
        基础信息
        '''
        return Text.assemble(
            'Developer: ', 
            ('YuHang Wu', 'bold'), 
            ' ( School of Pharmaceutical Sciences, Xiamen University )',
            '\nEmail: ',('yuhangxmu@stu.xmu.edu.cn', 'purple'),
            '\nGithub: https://github.com/CyberCatQ/pyCADD'
            )

    @property
    def version_info(self) -> None:
        '''
        版本信息
        '''
        return Text.assemble(
            'Version:  ',
            (self.version, 'bold blue'), 
            ' ' * (9-len(self.version)),
            'Last update:  ', 
            (self.update_date, 'bold blue')
            )

    @property
    def system_info(self) -> None:
        '''
        系统基本信息
        '''
        return Text.assemble(
            'Platform: ', 
            (platform.system(), 'bold blue'),
            ' ' * (9-len(platform.system())),
            'Current date: ', 
            (now, 'bold blue'),
            '\nNumber of parallel threads: ', 
            (str(os.cpu_count()), 'bold blue'), 
            ' Python Version: ',(platform.python_version(), 'bold blue'),
            '\nSchrodinger: ', self.schrodinger,
            '\nGaussian: ', self.gauss_dir
            )

    @property
    def info_index(self):
        self._info_index += 1
        return self._info_index

    def create_panel(self, options: list = None, additional_info:'str | dict'=None, options_label: str = 'Analysis Options', show_panel: bool = True) -> None:
        '''
        建立并渲染UI
        Parameters
        ----------
        options : list
            选项框内容
        additional_info : str | dict
            选项框上方的额外信息  
            传入字典时 可用于修改已存在的同key信息内容
        options_label : str
            选项框标签名
        show_panel : bool
            是否显示UI
        '''
        if options:
            self.options = options
        else:
            options = self.options

        if isinstance(additional_info, str):
            self._additional_info_dict[str(self.info_index)] = additional_info
            # self.additional_info = self.additional_info + '\n' + additional_info
        elif isinstance(additional_info, dict):
            for key, value in additional_info.items():
                self._additional_info_dict[key] = value

        grid_upper = Table(Column(self.title, justify='center'),
                           expand=True, show_edge=False, box=box.SIMPLE, padding=(1, 1))

        #grid_upper.add_row(Padding(self.title, 1))
        #grid_upper.add_row('-' * 48)

        grid_mid = Table.grid(expand=True)
        grid_mid.add_column(justify='center')
        grid_mid.add_row(self.version_info)
        grid_mid.add_row(self.system_info)
        grid_mid.add_row(self.basic_info)

        grid_lower = Table.grid(expand=True, padding=(0, 3, 0, 3))
        grid_lower.add_column(justify='left')
        grid_lower.add_column(justify='left')

        if options:
            left_num = len(options) // 2
            right_num = len(options) - left_num
            for i in range(right_num):
                try:
                    left_index = i
                    right_index = i + right_num
                    left = options[left_index]
                    right = options[right_index]
                except IndexError:
                    right = ''
                grid_lower.add_row(left, right)

        if self._additional_info_dict:
            info_text = ''.join(info + '\n' for info in self._additional_info_dict.values()).strip()
            additional_column = Padding(
                '[bold]' + info_text, (1, 0, 0, 3))
        else:
            additional_column = ''

        self.panel = Panel(Group(
            grid_upper,
            grid_mid,
            additional_column,
            Padding(self.menu_name, (1, 0, 0, 3)),
            Panel(
                grid_lower, title='[bold]%s' % options_label, title_align='left', padding=(1, 2))), expand=False)
        
        if show_panel:
            print(self.panel)
    
    def get_input(self, text:str, choices:list=None, default=None, show_default:bool=True, show_choices:bool=False):
        '''
        读取输入指令 返回flag
        '''
        return Prompt.ask(text, choices=choices, default=default, show_choices=show_choices, show_default=show_default)
    
    def get_confirm(self, text:str, default=True):
        '''
        读取输入指令 返回确认值
        '''
        return Confirm.ask(text, default=default)
    
    def clear_info(self):
        '''
        清空额外信息内容
        '''
        self.additional_info = ''
