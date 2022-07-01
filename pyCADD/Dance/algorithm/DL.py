from pandas import DataFrame
import torch
from torch import nn


class MLP(nn.Module):
    '''
    pytorch实现的多层感知机
    '''

    def __init__(self, input_dim: int, hidden_dim: int, output_dim: int, device: torch.device = None) -> None:
        '''
        初始化

        Parameters
        ----------
        input_dim : int
            输入维度
        hidden_dim : list
            隐藏层维度
        output_dim : int
            输出维度
        device : torch.device
            模型使用的设备
        '''
        super(MLP, self).__init__()
        self.input_dim = input_dim
        self.hidden_dim = hidden_dim
        self.output_dim = output_dim
        self.device = device if device is not None else torch.device(
            'cuda' if torch.cuda.is_available() else 'cpu')

        # 隐藏层数
        self.layer_num = len(hidden_dim)

        # 定义模型
        model = nn.Sequential()
        model.add_module('input_layer', nn.Linear(input_dim, hidden_dim[0]))
        model.add_module('relu_1', nn.ReLU())

        for layer_idx in range(self.layer_num - 1):
            model.add_module(
                f'hidden_layer_{layer_idx + 1}',
                nn.Linear(hidden_dim[layer_idx], hidden_dim[layer_idx + 1])
            )
            model.add_module(
                f'relu_{layer_idx + 2}',
                nn.ReLU()
            )

        model.add_module(
            'output_layer',
            nn.Linear(hidden_dim[-1], output_dim)
        )

        self.model = model.to(device)

    def forward(self, x):
        return self.model(x)

    def get_params(self):
        return {
            'device': self.device,
            'input_dim': self.input_dim,
            'output_dim': self.output_dim,
            'hidden_dim': self.hidden_dim,
        }

    def initialize(self, method='normal'):
        '''
        初始化模型参数

        Parameters
        ----------
        method : str
            初始化方法
            normal : 正态分布初始化
            xavier_uniform : Xavier 均匀初始化
            xavier_normal : Xavier 正态分布初始化
            kaiming_uniform : Kaiming 均匀初始化
            kaiming_normal : Kaiming 正态分布初始化
        '''
        if method == 'xavier_uniform':
            init_func = nn.init.xavier_uniform_
        elif method == 'xavier_normal':
            init_func = nn.init.xavier_normal_
        elif method == 'normal':
            init_func = nn.init.normal_
        elif method == 'kaiming_uniform':
            init_func = nn.init.kaiming_uniform_
        elif method == 'kaiming_normal':
            init_func = nn.init.kaiming_normal_
        else:
            raise ValueError('Unknown initialization method')
        for layer in self.model:
            if isinstance(layer, nn.Linear):
                init_func(layer.weight)

    def predict(self, X: DataFrame):
        '''
        推断

        Parameters
        ----------
        X : DataFrame
            输入数据

        Returns
        -------
        numpy.array
            推断结果
        '''
        self.eval()
        X = torch.tensor(X.values, dtype=torch.float)
        with torch.no_grad():
            X = X.to(self.device)
            y = self.forward(X)
            y = y.cpu().argmax(dim=1).numpy()
        return y

    def predict_proba(self, X: DataFrame):
        '''
        推断概率

        Parameters
        ----------
        X : DataFrame
            输入数据

        Returns
        -------
        numpy.array
            各分类推断概率
        '''
        self.eval()
        X = torch.tensor(X.values, dtype=torch.float)
        with torch.no_grad():
            X = X.to(self.device)
            y = self.forward(X)
            y = y.cpu().softmax(dim=1).numpy()
        return y
