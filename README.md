pyCADD
==========

## Features

* 调用 Gaussian 计算、分析配体分子结构优化、单点能、RESP(2)电荷
* 自动运行 AMBER 分子动力学模拟准备、运行和 MD 轨迹的基本分析
* 自动化调用Schrodinger Python API执行晶体准备、格点文件生成与对接、MMGBSA结合能计算等功能
* 调用多核并行计算 实现集合式对接与结果提取、数据分析
* CLI 快速调用

## Platform  

* Linux  

## Requirements 

### Python Version  

* 3.10 or Higher  

### Environment Requirements

依赖库详见 [environment.yml](environment.yml) 文件。
* 建议使用 Conda/Mamba 创建独立的Python环境以避免包冲突

```bash
conda create -n pyCADD python=3.10
conda activate pyCADD
conda install -f environment.yml
pip install pyCADD
```

### Additional Requirements

pyCADD的不同模块需要安装不同的软件来完成自动化工作流。如您仅需要使用其中的部分功能, 可根据您的需求仅安装相应的软件。  
当模块缺少必要软件时，您可能会在控制台看到相应的提示信息。

| Module   |  Software   | Version |
| -----------------  | -------- | -------- |
| Dock     | [Schrodinger Suite](https://www.schrodinger.com/) | 2020-3 or newer |
| Dynamic  | [AMBER](http://ambermd.org/) | 22 or newer |
| Density & Dynamic | [Gaussian](http://gaussian.com/) | 16.A01 or newer(optional) |
| Dynamic | [CUDA](https://developer.nvidia.com/cuda-zone) | 11.7 or newer (optional for `pmemd.cuda`) |

### ！Attention
* `pyCADD` 不包含以上所需程序的安装与许可证 您需要自行获得授权并安装恰当
* 使用本应用程序进行学术研究必须遵守以上所需各程序的相关文献引用规定

## Installation

`pyCADD`已发布至PyPI, 使用命令

    conda install -f environment.yml
    pip install pyCADD

即可安装pyCADD。

为了便于从命令行快速调用pyCADD的模块，提供了以下CLI接口：

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

## Usage
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
此工具包仅限于学习和批评使用, 维护周期不定，不建议用于商业用途。 

YH. W  
School of Pharmaceutical Sciences, Xiamen University  
Copyright © 2025 XMU   
2025-10-28

## 更新日志
* 2.0.0 (2025-10-28)
  * 重构 Density 模块
    * 无需依赖 Multiwfn 即可完成 RESP2 计算结果的解析与分析
    * 明确模块功能定位，去除不常用功能
    * 完全重写模板文件生成逻辑，去除第三方代码依赖
  * 重构 Dynamic 模块，采用更加模块化的设计，提升代码可维护性
    * 不再限制配体分子结构输入格式为 pdb
  * 重构 Dock 模块，提升代码可维护性
  * 所有模块 UI 功能停止维护并将逐渐移除，改为主要基于 CLI 调用
  * 修复了一些已知 BUG
  * 为所有关键模块重写了英文标准化代码注释
  * 更新了部分文档
  * 安装配置文件更新为 pyproject.toml，移除 setup.py

* 1.6.9 (2023-12-26)
  * 将 Demand 数据获取修改为直接与PDB服务器API交互，确保获取最新PDB数据
  * 添加了对PDB蛋白链突变信息的解析
  * 为 Dock 添加了多线程批量下载PDB的接口 `pycadd-dock download`
  * 其他 BUG 修复

* 1.6.8 (2023-11-21)
  * 为 Dock 模块 CLI 接口 pycadd-dock 添加了自动创建 schrodinger 虚拟环境的功能
  * 为 Dynamic 模块增加了额外的必须软件检查
  * 调整了MD分析工具中能量分解参数的传递方式，现在起始、结束、步长参数通过命令行直接传递
  * 增加了 mpirun 命令的 localhost:N 参数以避免在某些系统中出现核心数量不匹配的问题

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
  * 添加了Dynamic模块的MD后处理工具Analyzer
  * 为MD后处理添加了CLI接口 pycadd-mdanalysis
  * 为Dance添加了新的统计方法 标准富集因子NEF

* 1.6.0 (2022-07-07)
  * 分子动力学模块 Dynamic 开发完成
  * 添加了Dynamic的CLI接口
  * 为项目文档增加了 Dynamic 的 User Guide

* 1.5.4 (2022-07-04)
  * 基于在 jupyter notebook 中的应用需要，重构了Dance模块，简化了相关逻辑
  * 移除了Dance的UI, 完全使用库导入的方式
  * 为项目文档增加了Dance的 Use Guide

* 1.5.3 (2022-06-12)
  * 修复了Dance模块的一些已知BUG

* 1.5.2(2022-06-01)
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