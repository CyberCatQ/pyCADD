pyCADD VSW模块
==============

## 功能特性
* 自动调用Schrodinger VSW API执行虚拟筛选工作流
* 通过配置文件管理和注册基因、化合物数据库 便于重复筛选
* 用户界面友好
  
## Platform  

* Linux  

## Required

* [Schrodinger Suite](https://www.schrodinger.com/)2018-2 或更高版本
  
### Python Version  

* 3.6 or Higher  

### Python Modules  

* rich
* pandas
  
如果缺少python库 请使用`pip install <package>`命令安装

## How to Use
使用 `git clone` 命令 将 `pyCADD` 文件夹移入python的site-packages目录下 (setuptools开发中)  
然后, 您可以使用 `python -m pyCADD` 在任意您想要保存项目文件的目录下执行应用程序  
VSW过程文件将会全部存储于`vsw/`目录中

### 基因注册  

首先准备好基因名称同名的PDB列表文件  
以RXRa为例  *RXRa.txt*

    3OAP,9CR
    4K6I,9RA

然后 请您使用该模块注册需要进行VSW的基因名称  

* 注册需要您输入基因名称和上述txt文件的路径。

### 化合物库注册  

首先准备好您的化合物库结构文件( **.mae* **.maegz*)
然后 请您使用该模块注册需要进行VSW的化合物库  

* 注册需要您输入自定义的化合物库名和上述结构文件的路径。

当您完成基因和化合物注册后 即在基因和化合物库选择中可见 并无需再在多次VSW中重复输入相应文件的路径
  
* * *
此脚本仅限于学习和批评使用, 请勿用作其他用途。  
源码仅包含中文注释。

YH. W  
School of Pharmaceutical Sciences, Xiamen University  
2021-12-15
