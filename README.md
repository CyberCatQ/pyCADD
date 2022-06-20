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
* 1.5.3 (2022-06-12)
  * 修复了Dance模块的一些已知BUG

* 1.5.2 (2022-06-01)
  * 修复MLP模型的BUG
  * 修复获取预测结构函数的BUG
  * 添加获取预测结构的CLI支持
  * 修复旧版本ConcurrentlogHandler导致的死锁BUG, 不再限制setuptools的版本

* 1.5.1 (2022-05-18) 
  * 完全重构query模块 输入方式修改为Uniprot ID
  * 增加CLI快速调用接口

* 1.5.0 (2022-05-12) 
  * 重新以面向对象化方式完全重构 Dock 模块 降低耦合度
  * 修改了输入文件逻辑 
  * 将 Multidock 与 Dock 合并为单一模块 移除了 Multidock 的用户界面(UI)
  * 增加了 Dock/Multidock 的命令行接口(CLI) 比起UI更加快速、便捷、易用
  * 增加了 Dance 模块的机器学习与深度学习模型支持及简易性能评估(基于AUC)
  * 增加了 Dance 模块的标准化数据集性能快速评估工作流(evaluate_workflow)

## Platform  

* Linux  

## Required

* [Schrodinger Suite](https://www.schrodinger.com/)2020-3 或更高版本
* [AMBER](http://ambermd.org/) 18 或更高版本
* [Gaussian](http://gaussian.com/) 16.A01 或更高版本
* [Multiwfn](http://sobereva.com/multiwfn/) 3.7 或更高版本
* CUDA 9.0 或以上版本
### ！Attention
* `pyCADD` 不包含以上所需程序的安装与许可证 您需要自行获得授权并安装恰当
* 使用本应用程序进行学术研究必须遵守以上所需各程序的相关文献引用规定
  
### Python Version  

* 3.8 or Higher  

<a id='python-modules'></a>
### Python Modules  

* pandas
* numpy
* rich
* ConcurrentLogHandler
* scikit-learn
* seaborn
* xlsxwriter  

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


| 开发中 | Function |  
|-------|-------|
|*pyMDprepare.py*    | 准备AMBER MD必要文件(拓扑及坐标、力场参数等文件)
|*RESP2.sh*          | 调用Gaussian执行坐标优化并计算RESP2(0.5)电荷(已集成于 *pyMDprepare.py*)  |
|*pyPMEMD.py*        | 调用AMBER PMEMD(GPU加速)执行能量最小化、体系加热与分子动力学模拟  |
|*pynalysis.py*      | MD轨迹分析 输出RMSD/RMSF、氢键、二面角等变化情况 提取最低势能构象  调用MMPBSA计算吉布斯自由能变/熵变  |

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

使用

    pycadd-dock --help

来获取更多帮助信息。

### **注意**  

* 由于Schrodinger Suite使用了闭源环境提供Python API调用, 使用 `Dock` 模块相关功能(包括CLI命令)脚本时, 应该保证您正处于由 Schrodinger 建立的虚拟环境中，并在该环境中安装包括 `pyCADD` 在内所有需要的[Python Modules](#python-modules)
* 要使用Schrodinger内部脚本创建虚拟环境, 请参阅 [如何建立Schrodinger Python Virtual Env](https://content.schrodinger.com/Docs/r2022-1/python_api/intro.html#general-python-information)
* 在Schrodinger内部环境中安装额外的python packages, 请使用`run python3 -m pip install <package>`命令或参阅[Schrodinger Python API 文档](https://www.schrodinger.com/pythonapi/)

更多帮助信息, 请参阅[API文档](https://cybercatq.github.io/pyCADD)。

## Quick Start

### Docking of Many-to-One Function

`pyCADD` UI 的 `1. Dock Mode` 提供了将一个或多个配体对接于单个受体结构的功能  
运行 `pyCADD` 选择模式1并使用功能6 提供包含一个多个配体的单一 `.mae` 或 `.maegz` 的文件路径即可

### Docking of One-to-Many Function

`pyCADD` 的 `pycadd-dock ensemble-dock` 接口提供了高性能多核并行集合式对接计算及数据处理  
请确保当前工作目录 (Current Working Directory) 在您想要保存项目文件的目录中, 并在目录中额外准备：

* 一个分行列出的, 受体蛋白所属PDB ID的列表文本文件 **.txt* 或 **.csv*
* 需要对接的外源配体 **单一3D结构文件** **.pdb*  **.mae* 或其他Schrodinger支持的格式 (可选)
  
以下是一个示例 *example.txt*：  
！如果晶体包含有多个不同名称的小分子配体, 请在文件中指明 (以逗号分隔)以定义对接格点的中心

    3OAP,9CR
    5JI0,BRL
    4K6I,9RA

以及一个准备就绪的化合物库文件(其中含有若干个结构) library_file.mae

然后 您可以通过命令
    
    pycadd-dock ensemble-dock example.txt library_file.mae

启动 Ensemble Docking。

使用 

    pycadd-dock ensemble-dock --help

以获取更多帮助信息。


### 注意：

* 脚本将自动识别文本文件中每一行的PDB ID并自行批量下载  
* 应用程序默认使用系统最大核心数量的 75% 进行集合式对接工作 通过设置 `-n/--parallel <NUM>` 参数来指定使用的核心数量
* 如您需要执行晶体共结晶配体的回顾性对接(Self-Docking), 请添加参数 `--redock`
* 所有集合式对接工作完成后, 将自动提取重要的对接结果数据, 并保存在 `result` 目录下的 `_DOCK_FINAL_RESULTS.csv`字样的文件中, 且将产生汇总矩阵文件`matrix.csv`

## Gaussian Calculation Module
您可以直接运行 `pycadd-gauss` 并依据提示输入初始结构文件路径

或使用命令行  
    
    pycadd-gauss [input_file_path]

来直接载入初始结构文件

* * *
`pyCADD` 完全基于开发者自身实际需求开发, 更多其他功能及模块的使用教程尚未撰写, 请参阅[API文档](https://cybercatq.github.io/pyCADD)。  
如果您有任何问题或建议, 请在[项目主页](https://github.com/CyberCatQ/pyCADD)提交Issue。  

此脚本仅限于学习和批评使用, 请勿用作其他用途。  
源码仅包含中文注释。

YH. W  
School of Pharmaceutical Sciences, Xiamen University  
2022-05-18
