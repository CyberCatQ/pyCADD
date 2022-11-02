pyCADD
==========

## 功能特性

* 自动化调用Schrodinger Python API执行晶体准备、格点文件生成与对接、MMGBSA结合能计算等功能
* 调用多核并行计算 实现集合式对接与结果提取、数据分析
* 调用Gaussian、Multiwfn计算、分析配体分子结构优化、单点能、RESP(2)电荷
* AMBER分子动力学模拟准备、运行和MD轨迹的基本分析
* 用户界面友好
* 支持CLI快速调用

## 更新日志
* 1.6.3 (2022-10-27)
  * 修改了 Dock 模块的对接函数调用 不再使用Schrodinger任务管理层
    * 显著降低无意义时间损耗 提升 Ensemble Docking 速度
    * 完全减除任务管理层perl的长期内存占用 避免后期速度下降
    * 完全保障了 CPU 100%时间满载运行对接任务
  * 修改 Dock 的默认的结果文件生成为仅输出配体结构 
    * 显著减少对接结果文件体积
  * 修复了Dynamic Analysis模块的BUG
    * 现在以重原子间距离计算氢键键长
  
* 1.6.2 (2022-07-21)
  * 修复了 Dynamic-Analyzer 的一些BUG
  * 修复了 python wheel build 的一些BUG
  * 更新了API文档

* 1.6.1 (2022-07-14)
  * 添加了 Dynamic 模块的MD后处理工具 Analyzer
  * 为 Dynamic-Analyzer 添加了CLI接口 pycadd-mdanalysis
  * 为 Dance 添加了新的统计方法 标准富集因子NEF
  
* 1.6.0 (2022-07-07)
  * 分子动力学模块 Dynamic 开发完成
  * 添加了Dynamic的CLI接口
  * 为项目文档增加了 Dynamic 的 User Guide

## Platform  

* Linux  

## Required

* [Schrodinger Suite](https://www.schrodinger.com/)2020-3 或更高版本
* [AMBER](http://ambermd.org/) 18 或更高版本 
* [Gaussian](http://gaussian.com/) 16.A01 或更高版本
* [Multiwfn](http://sobereva.com/multiwfn/) 3.7 或更高版本
* [OpenBabel](https://openbabel.org/) 2.4 或更高版本
* CUDA 9.0 或以上版本(Optional)
### ！Attention
* `pyCADD` 不包含以上所需程序的安装与许可证 您需要自行获得授权并安装恰当
* 使用本应用程序进行学术研究必须遵守以上所需各程序的相关文献引用规定
  
### Python Version  

* 3.9 or Higher  

<a id='python-modules'></a>
### Python Modules  

* pandas
* numpy
* rich
* concurrent_log_handler
* scikit-learn
* seaborn
* xlsxwriter  
* pyyaml
* click
* scipy

#### NOTE:

使用`pip install pyCADD`命令时 将自动安装所需的packages

您也可以使用`pip install <package>`命令安装

## pyCADD Function

|        Name        | Function |
| -----------------  | -------- |
|*Dock* | 自动执行PDB晶体获取、优化、格点文件生成、对接等命令 |
|*Dance*    | 数据处理脚本集合  |  
|*VSW*            | 自动化虚拟筛选 |
|*query* | 自动化晶体信息查询与解析 |
|*Gauss* | 自动编写高斯结构优化、单点能、TDDFT等计算任务输入文件并启动计算任务 |
|*Dynamic* | 自动化分子动力学模拟准备与运行 自动化分子动力学模拟后处理分析 |

## How to Use

`pyCADD`已发布至PyPI

使用命令

    pip install pyCADD

即可安装 `pyCADD`  

随后 您可以使用命令 `pycadd` 或 `pyCADD` 来启动应用程序  
`pyCADD` 提供一个用户友好的界面使您能够轻易使用需要的功能, 请自行尝试。

为了便于从命令行快速调用pyCADD的模块，还提供了以下CLI接口：

| CLI | Module |
| ---- | ------ |
| `pycadd-dock` | *Dock* |
| `pycadd-query` | *query* |
| `pycadd-gauss` | *Gauss* |
| `pycadd-dynamic` | *Dynamic* |
| `pycadd-dynamic analysis` | *Dynamic-Analyzer* |

使用如

    pycadd-dock --help
    pycadd-dynamic --help
    pycadd-dock ensemble-dock --help

样式的命令来获取各模块的更多帮助信息。

各模块使用指南及更多帮助信息, 请参阅[API文档](https://cybercatq.github.io/pyCADD)。

* * *
`pyCADD` 完全基于开发者自身实际需求开发, 更多其他功能及模块的使用教程尚未撰写, 请参阅[API文档](https://cybercatq.github.io/pyCADD)。  
如果您有任何问题或建议, 请在[项目主页](https://github.com/CyberCatQ/pyCADD)提交Issue。  

此脚本仅限于学习和批评使用, 请勿用作其他用途。  
源码仅包含中文注释。

YH. W  
School of Pharmaceutical Sciences, Xiamen University  
Copyight © 2022 XMU   
2022-10-27
