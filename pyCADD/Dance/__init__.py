'''
Data Analyzer for Computer-aid drug design
'''

import json
import os
import pickle
from time import sleep

import hiddenlayer as hl
import numpy as np
import torch
from pandas import DataFrame
from pyCADD.Dance.algorithm import (Average, Consensus, Geo_Average, Minimum,
                                    MyMLP)
from pyCADD.Dance.core import (CV_model_evaluation, hyperparam_tuning,
                               split_data)
from sklearn.dummy import DummyClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (accuracy_score, classification_report,
                             confusion_matrix, f1_score, precision_score,
                             recall_score, roc_auc_score)
from sklearn.model_selection import train_test_split
from sklearn.naive_bayes import GaussianNB
from torch import nn
from torch.utils.data import DataLoader, TensorDataset
from torch.utils.data.sampler import WeightedRandomSampler
from xgboost import XGBClassifier


class Evaluator:

    def __init__(self, data:DataFrame, test_size:float=0.2, random_state:int=0):
        '''
        Parameters
        ----------
        data : DataFrame
            数据集
        test_size : float
            测试集比例
        random_state : int
            随机种子
        '''
        self.data = data
        self.test_size = test_size
        self.random_state = random_state
        self.log_step_interval = 10

        self.cv_results = None
        self.best_model = None

        self.classifiers = {}
        self.gbt_estimator_hyperparameters = {
        'gpu_id' : [0],
        'tree_method' : ['gpu_hist'],
        'n_estimators' : [200, 250, 300],
        'max_depth'     : [3, 5, 10, 20],
        'gamma'        : [0.01, 0.1, 0.5, 1],
        'learning_rate': [0.01, 0.05, 0.1],
        'subsample'    : [0.3, 0.5, 0.6],
        'alpha'        : [0.01, 0.1, 0.5, 1],
        'colsample_bytree': [0.3, 0.5, 1],
        'objective'    : ['binary:logistic'],
        'eval_metric': ['auc'],
        'use_label_encoder': [False]
        }
        
        self.lr_estimator_hyperparameters = {
        'C' : [0.01, 0.1, 1],
        'penalty': ['l1', 'l2']
        }   

        self.rf_estimator_hyperparameters = {
        'n_estimators' : [200, 250, 300],
        'max_depth'     : [3, 5, 10, 20],
        'max_features'  : ['auto', 'sqrt', 'log2'],
        'min_samples_split' : [2, 5, 10],
        'min_samples_leaf' : [1, 2, 5],
        }

        self._preprocess_data()

    @property
    def cs_classifiers(self):
        '''
        通识算法分类器
        '''
        return {
        'cs_mean': Average(lower_is_better=True),
        'cs_GeometricMean': Geo_Average(lower_is_better=True),
        'cs_Min': Minimum(lower_is_better=True)
        }

    @property
    def mlp_parameters(self):
        '''
        多层感知机参数 所有参数均在此调节
        '''
        return {
            'batch_size' : 128,
            'epochs' : 1000,
            'device': torch.device('cuda:0' if torch.cuda.is_available() else 'cpu'),
            'input_dim': self.X_train.shape[1],
            'output_dim' : len(set(self.y_train)),
            'hidden_dim' : (32, 16),
            'lr' : 0.01,
            'weight_decay' : 1e-4,
            'dropout1' : 0.2,
            'dropout2' : 0.5,
            'min_loss' : 0.1,
            'min_epochs': 100
        }

    @staticmethod
    def get_weights(y):
        '''
        计算不平衡数据集中各类标签的权重
        '''
        # 样本不平衡时，采用权重采样
        label_to_count = y.value_counts()
        _weights = []
        for label in y:
            _weights.append(1 / label_to_count[label])
        return _weights

    @property
    def train_dataset(self):
        '''
        训练集
        '''
        return TensorDataset(torch.tensor(self.X_train.values, dtype=torch.float), torch.tensor(self.y_train.values, dtype=torch.long))
    
    @property
    def train_dataloader(self, weighted:bool=True):
        '''
        训练集数据加载器 按照权重采样
        '''
        weights = self.get_weights(self.y_train)
        sampler = WeightedRandomSampler(weights, len(weights)) if weighted else None
        batch_size = self.mlp_parameters['batch_size']
        return DataLoader(self.train_dataset, batch_size=batch_size, sampler=sampler)
    
    @property
    def test_dataset(self):
        '''
        测试集
        '''
        return TensorDataset(torch.tensor(self.X_test.values, dtype=torch.float), torch.tensor(self.y_test.values, dtype=torch.long))
    
    @property
    def test_dataloader(self, weighted:bool=True):
        '''
        测试集数据加载器 按照权重采样
        '''
        weights = self.get_weights(self.y_test)
        sampler = WeightedRandomSampler(weights, len(weights)) if weighted else None
        batch_size = self.mlp_parameters['batch_size']
        return DataLoader(self.test_dataset, batch_size=batch_size, sampler=sampler)
    
    @property
    def mlp_classifier(self):
        '''
        多层感知机MLP分类器 参数于mlp_parameters调节
        '''
        parameters = self.mlp_parameters
        input_dim = parameters['input_dim']
        output_dim = parameters['output_dim']
        hidden_dim = parameters['hidden_dim']
        dropout1 = parameters['dropout1']
        dropout2 = parameters['dropout2']
        lr = parameters['lr']
        weight_decay = parameters['weight_decay']
        batch_size = parameters['batch_size']
        epochs = parameters['epochs']

        return MyMLP(input_dim, hidden_dim, output_dim, dropout1, dropout2, lr, weight_decay, batch_size, epochs)

    # 定义处理数据的标准工作流程
    def _preprocess_data(self) -> None:
        '''
        划分数据集为特征和标签 并以0填充缺失值
        '''
        X, y = split_data(self.data, 'activity', True, 1)
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=self.test_size, random_state=self.random_state)
        
        self.X_train = X_train
        self.y_train = y_train
        self.X_test = X_test
        self.y_test = y_test

        return X_train, X_test, y_train, y_test
    
    def get_parameters(self, classifier_name, estimator=None, hyperparameters=None, method='random', force_search=False) -> dict:
        '''
        加载或搜索分类器参数  

        存在best_params_[classifier_name].json时, 加载参数 否则搜索参数

        Parameters
        ----------
        classifier_name : str
            分类器名称
            lr : LogisticRegression
            rf : RandomForestClassifier
            gbt : XGBClassifier
        
        estimator : sklearn.base.BaseEstimator
            分类器实例
        hyperparameters : dict
            分类器参数搜索范围
        method : str
            搜索方法['random', 'grid']
        force_search : bool
            是否强制重新进行网格搜索

        Returns
        -------
        dict
            分类器参数
        '''
        if estimator is None:
            if classifier_name == 'lr':
                clf_name = 'LogisticRegression'
                estimator = LogisticRegression(max_iter=1000)
                hyperparameters = self.lr_estimator_hyperparameters
            elif classifier_name == 'rf':
                clf_name = 'RandomForestClassifier'
                estimator = RandomForestClassifier()
                hyperparameters = self.rf_estimator_hyperparameters
            elif classifier_name == 'gbt':
                clf_name = 'XGBClassifier'
                estimator = XGBClassifier()
                hyperparameters = self.gbt_estimator_hyperparameters
            else:
                raise ValueError('Estimator is not specified.')
        else:
            if hyperparameters is None:
                raise ValueError('Hyperparameter grid is not specified.')
            clf_name = estimator.__class__.__name__

        # 已存在参数且不强制搜索 则直接读取
        if os.path.exists(f'best_params_{clf_name}.json') and not force_search:
            with open(f'best_params_{clf_name}.json', 'r') as f:
                parameters = json.load(f)
        # 否则进行超参搜索
        else:
            print(f'{classifier_name} parameters are not found, start searching...')
            parameters = hyperparam_tuning(estimator, hyperparameters, self.X_train, self.y_train, save=True, model_name=clf_name, method=method)
        return parameters

    def add_classifier(self, classifier_name, parameters=None, classifier_=None) -> None:
        '''
        添加分类器到分类器列表

        classifier_name : str
            预设分类器名称
            lr : LogisticRegression
            rf : RandomForestClassifier
            gbt : XGBClassifier
            dummy : DummyClassifier
            nbc : GaussianNB
        parameters : dict | None
            分类器参数 将会被加载到分类器实例中
        classifier_ : sklearn.base.BaseEstimator
            其他非预设的分类器实例
        '''
        if classifier_name == 'lr':
            classifier = LogisticRegression(max_iter=1000)
            clf_name = 'LogisticRegression'
        elif classifier_name == 'rf':
            classifier = RandomForestClassifier()
            clf_name = 'RandomForestClassifier'
        elif classifier_name == 'gbt':
            classifier = XGBClassifier()
            clf_name = 'XGBClassifier'
        elif classifier_name == 'dummy':
            classifier = DummyClassifier(strategy='stratified', random_state=42)
            clf_name = 'DummyClassifier'
        elif classifier_name == 'nbc':
            classifier = GaussianNB()
            clf_name = 'NaiveBayesClassifier'
        elif classifier_name == 'cs':
            self.classifiers.update(self.cs_classifiers)
            return None
        elif classifier_name == 'mlp':
            classifier = self.mlp_classifier
            clf_name = 'MLPClassifier'
        else:
            if classifier_ is None:
                raise ValueError('Classifier is required')
            else:
                classifier = classifier_
                clf_name = classifier_name.upper()

        if parameters is not None:
            try:
                classifier.set_params(**parameters)
            except AttributeError:
                pass
        self.classifiers[clf_name] = classifier
    
    def cross_validation(self, plot=True, **kwargs):
        '''
        将分类器列表中的分类器执行30x4交叉验证

        Parameters
        ----------
        plot : bool
            是否绘制交叉验证结果图表
        kwargs : dict
            其他传递给CV_model_evaluation的参数
        '''
        results = CV_model_evaluation(self.classifiers, self.X_train, self.y_train, plot=plot, **kwargs)
        report_results = {}
        for clf_name, clf in self.classifiers.items():
            report_results[clf_name] = np.mean(results[clf_name])
        
        report_results['Best_SCP'] = np.max(results['SCP'])
        report_results['Mean_SCP'] = np.mean(results['SCP'])
        report_results['Worst_SCP'] = np.min(results['SCP'])
        report_results['SCP_std'] = np.std(results['SCP'])

        return report_results

    # 定义报告函数
    @staticmethod
    def evaluate(y_test, y_pred, y_proba, print_=True):
        '''
        评估分类器
        Parameters
        ----------
        y_test : numpy.ndarray
            测试集标签
        y_pred : numpy.ndarray
            预测集标签
        y_proba : numpy.ndarray
            预测集概率
        print_ : bool
            是否打印评估结果
        '''
        acc = accuracy_score(y_test, y_pred)
        precision = precision_score(y_test, y_pred)
        recall = recall_score(y_test, y_pred)
        f1 = f1_score(y_test, y_pred)
        auc = roc_auc_score(y_test, y_proba)

        if print_:
            print(classification_report(y_test, y_pred))
            print('Confusion Matrix\n', confusion_matrix(y_test, y_pred))
            print('Accuracy:', acc)
            print('Precision:', precision)
            print('Recall:', recall)
            print('F1 score:', f1)
            print('ROC AUC:', auc)
        
        return {
            'accuracy': acc,
            'precision': precision,
            'recall': recall,
            'f1': f1,
            'auc': auc
        }
    
    def report(self):
        '''
        打印所有模型的评估结果
        '''
        for clf_name, clf in self.classifiers.items():
            print('-'*50)
            print('Model:', clf_name, '\n')
            if isinstance(clf, nn.Module) or isinstance(clf, Consensus):
                print('Skip')
                continue
            clf.fit(self.X_train, self.y_train)
            y_pred = clf.predict(self.X_test)
            y_proba = clf.predict_proba(self.X_test)[:, 1]
            self.evaluate(self.y_test, y_pred, y_proba, print_=True)
            print('-'*50)
    
    @staticmethod
    def init_weights(model):
        '''
        初始化模型权重
        '''
        for param in model.parameters():
            nn.init.normal_(param, mean=0, std=0.01)

    @staticmethod
    def train_model(model, optimizer, loss_function, train_dataloader, test_dataloader, epochs=500, device=None, min_epochs=100, min_loss=0.1, plot=False):
        '''
        训练模型

        保存最佳模型(Test loss低于min_loss, 训练epoch数大于min_epochs)参数

        Parameters
        ----------
        model : torch.nn.Module
            模型
        optimizer : torch.optim.Optimizer
            优化器
        loss_function : torch.nn.modules.loss
            损失函数
        train_dataloader : torch.utils.data.DataLoader
            训练集数据加载器
        test_dataloader : torch.utils.data.DataLoader
            测试集数据加载器
        epochs : int
            训练次数
        device : torch.device
            训练设备 CPU/GPU
        min_epochs : int
            最小训练次数
        min_loss : float
            最小损失
        plot : bool
            是否在训练过程中绘制状态(Loss/Acc/Auc)曲线

        Returns
        -------
        tuple[history, best_checkpoint]

        history : Hiddenlayer.History
            训练状态记录器 包含keys: train_loss, train_acc, train_auc, test_loss, test_acc, test_auc
            通过history[key].data获取数据
        best_checkpoint : torch.nn.Module
            最佳模型检查点参数字典

            {
                'model': model.state_dict(), 
                'optimizer': optimizer.state_dict(),
                'epoch': epoch + 1
            }
        '''
        if device is None:
            device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
        history = hl.History()
        best_checkpoint = None
        for epoch in range(epochs):
            
            train_loss = []
            train_acc = []
            train_auc = []
            model.train()
            for i, (train_data_batch, train_label_batch) in enumerate(train_dataloader):
                train_data_batch_gpu = train_data_batch.float().to(device)
                train_label_batch_gpu = train_label_batch.to(device)
                train_outputs = model(train_data_batch_gpu)

                loss = loss_function(train_outputs, train_label_batch_gpu)
                optimizer.zero_grad()
                loss.backward()
                optimizer.step()

                train_loss.append(loss.detach().cpu())
                train_outputs = train_outputs.detach().cpu()
                train_pred = train_outputs.argmax(dim=1)
                train_proba = train_outputs.softmax(dim=1)[:, 1]

                train_acc.append(accuracy_score(train_label_batch.numpy(), train_pred.numpy()))
                train_auc.append(roc_auc_score(train_label_batch.numpy(), train_proba.numpy()))

            train_loss = np.mean(train_loss)
            train_acc = np.mean(train_acc)
            train_auc = np.mean(train_auc)

            model.eval()  
            with torch.no_grad():
                test_acc = []
                test_auc = []
                test_loss = []

                for i, (test_data_batch, test_label_batch) in enumerate(test_dataloader):
                    test_data_batch_gpu = test_data_batch.float().to(device)
                    test_outputs = model(test_data_batch_gpu)
                    test_outputs = test_outputs.detach().cpu()

                    batch_test_loss = loss_function(test_outputs, test_label_batch)
                    test_pred = test_outputs.argmax(dim=1)
                    test_proba = test_outputs.softmax(dim=1)[:,1]

                    test_acc.append(accuracy_score(test_label_batch.numpy(), test_pred.numpy()))
                    test_auc.append(roc_auc_score(test_label_batch.numpy(), test_proba.numpy()))
                    test_loss.append(batch_test_loss.item())

                test_acc = np.mean(test_acc)
                test_auc = np.mean(test_auc)
                test_loss = np.mean(test_loss)

                if test_loss < min_loss and epoch >= min_epochs:
                    min_loss = test_loss
                    best_checkpoint = {'model': model.state_dict(), 'optimizer': optimizer.state_dict(),'epoch': epoch + 1, 'test_loss': test_loss}
            
            history.log(
                    epoch + 1, 
                    train_loss=loss.item(),
                    test_loss=test_loss.item(), 
                    train_acc=train_acc, 
                    train_auc=train_auc, 
                    test_acc=test_acc, 
                    test_auc=test_auc
                    )
                    
            if plot:
                canvas = hl.Canvas()
                with canvas:
                    canvas.draw_plot([history['train_loss'], history['test_loss']], ylabel='Loss')
                    canvas.draw_plot([history['train_acc'], history['test_acc']], ylabel='Accuracy')
                    canvas.draw_plot([history['train_auc'], history['test_auc']], ylabel='AUC')
            else:
                print(f'Epoch {epoch + 1}/{epochs} | Train Loss: {train_loss:.4f} | Test Loss: {test_loss:.4f} | Train Acc: {train_acc:.4f} | Test Acc: {test_acc:.4f} | Train AUC: {train_auc:.4f} | Test AUC: {test_auc:.4f}')
        
        return history, best_checkpoint

    def mlp_train(self, plot=True):
        '''
        基于MLP模型训练

        训练完成后 history 和 best_checkpoint 属性会被赋值于evaluator实例history和best_checkpoint属性
        最终模型(完成所有epoch训练)检查点保存到evaluator实例final_model属性

        Parameters
        ----------
        plot : bool
            是否在训练过程中绘制状态(Loss/Acc/Auc)曲线
        
        Return
        ------
        tuple(nn.Module, nn.Module)
            最佳模型, 最终模型
        '''
        X_train, X_test, y_train, y_test = self.X_train, self.X_test, self.y_train, self.y_test
        device = self.mlp_parameters['device']
        epochs = self.mlp_parameters['epochs']
        weight_decay = self.mlp_parameters['weight_decay']
        lr = self.mlp_parameters['lr']
        min_epochs = self.mlp_parameters['min_epochs']
        min_loss = self.mlp_parameters['min_loss']

        # 实例化一个模型
        model = self.mlp_classifier
        model.to(device)
        model.apply(self.init_weights)

        train_dataloader = self.test_dataloader
        test_dataloader = self.test_dataloader

        loss_function = nn.CrossEntropyLoss()
        # NOTE: 务必在初始化模型完成后再定义优化器
        optimizer = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=weight_decay)

        with torch.no_grad():
            current_auc = roc_auc_score(y_train, np.argmax(model(torch.tensor(X_train.values, dtype=torch.float).to(device)).detach().cpu().numpy(), axis=1))
        print('Current AUC:', current_auc)
        print('Start Training...')
        sleep(1)
        
        history, best_checkpoint = self.train_model(model, optimizer, loss_function, train_dataloader, test_dataloader, epochs, device, min_epochs, min_loss, plot)
        self.history = history
        self.best_checkpoint = best_checkpoint
        self.final_checkpoint = {'model': model.state_dict(), 'optimizer': optimizer.state_dict(), 'epoch': epochs}

        self.best_model = self.mlp_classifier
        self.best_model.load_state_dict(best_checkpoint['model'])
        self.final_model = model

        return self.best_model, self.final_model

    def print_classifier_info(self):
        '''
        打印分类器参数信息
        '''
        formatter = '{0:<25}{1:<30}{2:<40}'
        print('='*100)
        print(formatter.format('Classifier', 'Parameters', 'Values'))
        print('-'*100)
        for clf_name, clf in self.classifiers.items():
            parameters = [(k,v) for k, v in clf.get_params().items()]
            print(formatter.format(clf_name, str(parameters[0][0]), str(parameters[0][1])))
            for k, v in parameters[1:]:
                print(formatter.format('', str(k), str(v)))
            print('-'*100)
        print('='*100)

    def print_cv_results(self):
        '''
        打印交叉验证结果
        '''
        formatter = '{0:<25}{1:<30}'
        print(formatter.format('Classifier', 'AUC'))
        print(formatter.format('-'*25, '-'*30))
        for k, v in self.cv_results.items():
            print(formatter.format(k, v))

    def draw_after_training(self):
        '''
        训练完成后 绘制训练过程状态曲线
        '''
        history = self.history
        canvas = hl.Canvas()
        with canvas:
            canvas.draw_plot([history['train_loss'], history['test_loss']], ylabel='Loss')
            canvas.draw_plot([history['train_acc'], history['test_acc']], ylabel='Accuracy')
            canvas.draw_plot([history['train_auc'], history['test_auc']], ylabel='AUC')

    def _save_model(self, model_state_dict, save_path):
        '''
        保存模型状态参数

        Parameters
        ----------
        model_state_dict : dict
            模型检查点状态字典
        save_path : str
            保存文件路径
        '''
        model_state_dict.update(
            {
                'structure': 
                {
                    'input_dim':self.mlp_parameters['input_dim'],
                    'hidden_dim': self.mlp_parameters['hidden_dim'],
                    'output_dim': self.mlp_parameters['output_dim'],
                    'dropout1': self.mlp_parameters['dropout1'],
                    'dropout2': self.mlp_parameters['dropout2'],
                }
            }
        )
        torch.save(model_state_dict, save_path)

    def save_model(self, dir_path='./models'):
        '''
        保存最佳模型和最终模型 并保存训练记录数据

        Parameters
        ----------
        dir_path : str
            保存目录路径
        '''
        self._save_model(self.best_checkpoint, os.path.join(dir_path, 'best_model.pth'))
        self._save_model(self.final_checkpoint, os.path.join(dir_path, 'final_model.pth'))
        with open(os.path.join(dir_path, 'history.pkl'), 'wb') as f:
            pickle.dump(self.history, f)

    def load_model(self, path):
        '''
        加载模型参数

        Parameters
        ----------
        path : str
            模型参数文件路径
        
        Returns
        -------
        tuple[nn.Module, torch.optim.Optimizer]
            加载状态后的MLP模型与优化器
        '''
        base_model = self.mlp_classifier
        checkpoint = torch.load(path)
        base_model.load_state_dict(checkpoint['model'])
        base_optimizer = torch.optim.Adam(params=base_model.parameters())
        base_optimizer.load_state_dict(checkpoint['optimizer'])

        return base_model, base_optimizer
    
    def load_history(self, path):
        '''
        加载训练记录数据 赋值于evaluator实例history属性
        '''
        with open(path, 'rb') as f:
            self.history = pickle.load(f)

    def evaluate_workflow(self, force_search=False):
        '''
        评估数据集的标准工作流
            具体执行以下步骤:

                1. 获取LR、XGBoost、RF三个分类器的超参数   
                2. 添加预设的分类器(LR、XGBoost、RF、NBC、Dummy)以及通识算法(Mean、GeometricMean、Minimum)到待评估列表   
                3. 训练MLP分类器   
                4. 对待评估列表中的分类器执行30x4交叉验证   
                5. 报告评估结果 并保存于 cv_results.json 文件中   

        '''
        lr_params = self.get_parameters('lr', method='grid', force_search=force_search)
        gbt_params = self.get_parameters('gbt', method='random', force_search=force_search)
        rf_params = self.get_parameters('rf', method='random', force_search=force_search)

        self.add_classifier('lr', parameters=lr_params)
        self.add_classifier('gbt', parameters=gbt_params)
        self.add_classifier('rf', parameters=rf_params)
        self.add_classifier('nbc')
        self.add_classifier('dummy')
        self.add_classifier('cs')
        self.add_classifier('mlp')
        self.print_classifier_info()

        print('='*100)
        print('Start Training MLP...')
        sleep(3)
        self.mlp_train(plot=False)
        self.save_model()
        print('Training finished. Use draw_after_training() to draw the results.')
        print('='*100)
        
        print('='*100)
        print('Start Cross-Validation...')
        sleep(1)
        self.cv_results = self.cross_validation()
        self.print_cv_results()
        with open('cv_results.json', 'w') as f:
            json.dump(self.cv_results, f)
        print('='*100)
