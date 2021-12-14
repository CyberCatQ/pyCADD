import os
from datetime import datetime

from rich import box
from rich.console import Group
from rich.padding import Padding
from rich.panel import Panel
from rich.table import Column, Table
from rich.text import Text

date = datetime.now()
year = str(date.year)
month = str(date.month)
day = str(date.day)
now = "%s-%s-%s" % (year, month, day)


class UI:
    '''
    pyCADD程序用户交互界面(user interface)
    '''

    def __init__(self, menu_name: str = 'Main', additional_info: str = '') -> None:
        self.version = '1.30'
        self.update_date = '2021-12-14'
        self.menu_name = '[bold magenta]Menu: %s' % menu_name
        self.additional_info = additional_info

    @property
    def title(self) -> None:
        '''
        程序标题样式
        '''
        return Text.assemble(
            ('pyCADD', 'bold medium_spring_green'), 
            ' -- A ', ('Python Script', 'bold yellow'), 
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
            ' ( School of Pharmaceutical Sciences, Xiamen University )'
            )

    @property
    def version_info(self) -> None:
        '''
        版本信息
        '''
        return Text.assemble(
            'Version ', 
            (self.version, 'bold blue'), 
            '   Last Update: ', 
            (self.update_date, 'bold blue')
            )

    @property
    def system_info(self) -> None:
        '''
        系统基本信息
        '''
        return Text.assemble(
            'Number of parallel threads:  ', 
            (str(os.cpu_count()), 'bold blue'), 
            '  Current date: ', 
            (now, 'bold blue')
            )

    def create_panel(self, options: list = None) -> None:
        '''
        布局定位
        '''

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

        if self.additional_info:
            additional_column = Padding(
                '[bold]' + self.additional_info, (1, 0, 0, 3))
        else:
            additional_column = ''

        self.panel = Panel(Group(
            grid_upper,
            grid_mid,
            additional_column,
            Padding(self.menu_name, (1, 0, 0, 3)),
            Panel(
                grid_lower, title='[bold]Analysis Options', title_align='left', padding=(1, 2))), expand=False)


class UI_dock(UI):
    '''
    单晶体对接UI
    '''

    def __init__(self, menu_name: str = 'Simple Mode', additional_info: str = '') -> None:
        super().__init__(menu_name=menu_name, additional_info=additional_info)


class UI_Multimode(UI):
    '''
    多晶体对接UI
    '''

    def __init__(self, menu_name: str = 'Multiple Mode', additional_info: str = '') -> None:
        super().__init__(menu_name=menu_name, additional_info=additional_info)
