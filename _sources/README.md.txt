pyCADD
==========

## 功能特性

* 自动化调用Schrodinger Python API执行晶体准备、格点文件生成与对接、MMGBSA结合能计算等功能
* 调用多核并行计算 实现集合式对接与结果提取、数据分析
* 调用Gaussian、Multiwfn计算、分析配体分子结构优化、单点能、RESP(2)电荷
* 自动运行 AMBER 分子动力学模拟准备、运行和 MD轨迹的基本分析
* 提供简单友好的用户界面
* 支持CLI快速调用

## 更新日志
* 1.6.7 (2023-11-15)
  * 修复了 Dynamic 解析输入文件过程的一个BUG，该BUG将导致部分限制性MD的输入文件解析失败
  * 调整了 Dynamic 在调用命令失败时的输出信息
  * 更新了部分文档和注释

* 1.6.6 (2023-11-10)
  * 重构 Dynamic 模块构建工作流的方式，现在更加模块化和自由
  * 修复了一些BUG
  * 更新了部分注释与文档
  * Dynamic 与 Density 导入时添加简单的必要检查

* 1.6.5 (2023-04-28)
  * 修复 Dock 开展 Ensemble Docking 完成时数据提取的BUG
  * 调整了一些代码的位置便于后续维护
  * 修复了一些其它 BUG
  
* 1.6.4 (2023-03-15)
  * Dynamic 支持自动化准备和模拟Apo蛋白
  * 为 Dynamic-Analyzer 添加了氢键lifetime分析功能
  * MD后处理分析默认工作流移除了提取最低能量构象的功能
  * 除非有log信息输出，否则 pyCADD 不再生成空log文件
  * Gauss 模块更名为 Density 模块
  * query 模块更名为 Demand 模块
  * 一些 Bug 修复
  
* 1.6.3 (2022-10-27)
  * 修改了 Dock 模块的对接函数调用 不再使用Schrodinger任务管理层
    * 显著降低无意义时间损耗 提升 Ensemble Docking 速度
    * 移除任务管理层perl的长期内存占用 避免后期速度下降
    * 保障了 CPU 100%时间满载运行对接任务
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
  * 为 Dynamic-Analyzer 添加了CLI接口 pycadd-dynamic analysis
  * 为 Dance 添加了新的统计方法 标准富集因子NEF

## Platform  

* Linux  

## Requirements 

### Python Version  

* 3.8 or Higher  

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

#### NOTE

使用 `pip install pyCADD` 命令时 将自动安装所需的packages

您也可以使用 `pip install <package>` 命令安装

### Additional Requirements

pyCADD的不同模块需要安装不同的软件来完成自动化工作流。如您仅需要使用其中的部分功能, 可根据您的需求仅安装相应的软件。  
当模块缺少必要软件时，您可能会在控制台看到相应的提示信息。

| Module   |  Software   | Version |
| -----------------  | -------- | -------- |
| Dock     | [Schrodinger Suite](https://www.schrodinger.com/) | 2020-3 or newer |
| Dynamic  | [AMBER](http://ambermd.org/) | 18 or newer |
| Dynamic | [OpenBabel](https://openbabel.org/)  | 2.4 or newer |
| Density & Dynamic | [Gaussian](http://gaussian.com/) | 16.A01 or newer |
| Dynamic | [Multiwfn](http://sobereva.com/multiwfn/) |3.7 or newer (optional for RESP charge) |
| Dynamic | [CUDA](https://developer.nvidia.com/cuda-zone) | 9.0 or newer (optional for pmemd.cuda) |

### ！Attention
* `pyCADD` 不包含以上所需程序的安装与许可证 您需要自行获得授权并安装恰当
* 使用本应用程序进行学术研究必须遵守以上所需各程序的相关文献引用规定

## Installation

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
| `pycadd-demand` | *Demand* |
| `pycadd-density` | *Density* |
| `pycadd-dynamic` | *Dynamic* |
| `pycadd-dynamic analysis` | *Dynamic-Analyzer* |

## pyCADD Function

|        Name        | Function |
| -----------------  | -------- |
|*Dock* | 自动执行PDB晶体获取、优化、格点文件生成、对接等命令 |
|*Dance*    | 数据处理脚本集合 |  
|*Demand* | 自动化晶体信息查询与解析 |
|*Density* | 自动编写高斯结构优化、单点能、TDDFT等计算任务输入文件并启动计算任务 |
|*Dynamic* | 自动化 AMBER MD模拟准备与运行 自动化MD模拟后处理分析 |

## How to Use
请参阅[文档](https://cybercatq.github.io/pyCADD)中的各模块使用指南。

或使用
```bash
pycadd-dock --help
pycadd-dynamic --help
pycadd-dock ensemble-dock --help
```
样式的命令来获取各模块的更多帮助信息。
各模块使用指南及更多帮助信息, 请参阅[API文档](https://cybercatq.github.io/pyCADD)。

* * *
`pyCADD` 完全基于开发者自身实际需求开发, 更多其他功能及模块的使用教程尚未完全撰写, 但可能附有详尽注释。可参阅[文档](https://cybercatq.github.io/pyCADD)API部分及代码注释。  
如果您有任何问题或建议, 请在[项目主页](https://github.com/CyberCatQ/pyCADD)提交Issue。  
如果任何模块对您有帮助，欢迎为此项目点亮星星。

UI等大部分功能为作者练手之作， 可能不具有实用性。  
此脚本仅限于学习和批评使用, 请勿用作其他用途。  
源码仅包含中文注释。

YH. W  
School of Pharmaceutical Sciences, Xiamen University  
Copyight © 2023 XMU   
2023-11-08
