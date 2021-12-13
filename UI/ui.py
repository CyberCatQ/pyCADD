from datetime import datetime

from rich import print
from rich.console import Console, Group
from rich.layout import Layout
from rich.padding import Padding
from rich.panel import Panel
from rich.table import Table, Column
from rich.text import Text

date = datetime.now()
year = str(date.year)
month = str(date.month)
day = str(date.day)


class UI:
    '''
    pyCADD程序用户交互界面(user interface)
    '''

    def __init__(self) -> None:
        pass

    @property
    def title(self) -> None:
        '''
        程序标题样式
        '''
        return Text.assemble(('Python Script', 'bold yellow'), ' For ', ('Computer-aid Drug Design', 'bold cyan'))

    def grid(self, options:list=None) -> None:
        '''
        布局定位
        '''
        grid_upper = Table.grid(expand=True)
        grid_upper.add_column(justify='center')
        grid_upper.add_row(self.title)
        grid_upper.add_row('-' * 48)
        grid_upper.add_row('Author: YH. W      Last Update: 2021-12-13')

        grid_lower = Table.grid(expand=True)
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
                
        self.window = Panel(Group(Panel(grid_upper), Panel(grid_lower)), expand=False)


class UI_dock(UI):
    '''
    单晶体对接UI
    '''

    def __init__(self) -> None:
        super().__init__()


class UI_Multimode(UI):
    '''
    多晶体对接UI
    '''
