pyCADD
==========

## 功能特性

* 自动化调用Schrodinger Python API执行晶体准备、格点文件生成与对接、MMGBSA结合能计算等功能
* 调用多核并行计算 实现集合式对接与结果提取、数据分析
* 调用Gaussian、Multiwfn计算、分析配体分子结构优化、单点能、RESP(2)电荷
* AMBER分子动力学模拟准备、运行和MD轨迹的基本分析
* 用户界面友好

## Platform  

* Linux  

## Required

* [Schrodinger Suite](https://www.schrodinger.com/)2018-2 或更高版本
* [AMBER](http://ambermd.org/) 18 或更高版本 (要使用PMEMD 请安装CUDA 9.0或以上版本)
* [Gaussian](http://gaussian.com/) 09.D01 或更高版本 (默认为g16)
* [Multiwfn 3.7](http://sobereva.com/multiwfn/)

### Python Version  

* 3.6 or Higher  

### Python Modules  

* rich
* pandas
* ConcurrentLogHandler
* pytraj  
* xlsxwriter  
* Openbabel

如果缺少python库 请使用`pip install <package>`命令安装

## pyCADD Function

|        Name        | Function |
| -----------------  | -------- |
|*Dock* | 自动执行PDB晶体获取、优化、格点文件生成、对接等命令(请使用$SCHRODINGER/run运行此脚本) |
|*Multidock* | 多核心并行计算集合式对接方案与数据提取 |
|*pyMDprepare.py*    | 准备AMBER MD必要文件(拓扑及坐标、力场参数等文件)
|*RESP2.sh*          | 调用Gaussian执行坐标优化并计算RESP2(0.5)电荷(已集成于 *pyMDprepare.py*)  |
|*pyPMEMD.py*        | 调用AMBER PMEMD(GPU加速)执行能量最小化、体系加热与分子动力学模拟  |
|*pynalysis.py*      | MD轨迹分析 输出RMSD/RMSF、氢键、二面角等变化情况 提取最低势能构象  调用MMPBSA计算吉布斯自由能变/熵变  |
|*dataprocess.py*    | 数据处理脚本集合  |  
|*vsw.py*            | 自动化虚拟筛选 |
|*libpynalysis/*        | pynalysis.py所需模块包|

## How to Use

使用 `git clone` 命令 将 `pyCADD` 文件夹移入python的site-packages目录下 (setuptools开发中)  
然后, 您可以使用 `python -m pyCADD` 在任意您想要保存项目文件的目录下执行应用程序  
`pyCADD` 提供一个用户友好的界面使您能够轻易使用需要的功能, 请自行尝试。

### **注意**  

* 通过 `python -m pyCADD` 使用本应用时 您无需额外操作
* 由于Schrodinger Suite使用了闭源环境提供Python API调用, 如需单独使用Schrodinger Suite相关功能脚本, 您需要确保运行环境处于Schrodinger内部环境中 并在该环境中安装需要的python packages
* 为了以Schrodinger内部环境运行任何应用程序 您需要在命令前添加 `run` 来进入该环境
* 在Schrodinger内部环境中安装额外的python packages, 请使用`run python3 -m pip install <package>`命令或参阅[Schrodinger Python API 文档](https://www.schrodinger.com/pythonapi/)

以PDBID 1XLS为例：  

    mkdir 1XLS
    cd 1XLS
    python -m pyCADD

    # 您也可以单独运行各模块
    # run python -m pyCADD.Dock
    # run python -m pyCADD.Multidock

更多帮助信息, 请参阅API文档。

## Schrodinger Multiple Docking

### Docking of Many-to-One Function

`pyCADD` 的 `Simple Mode` 提供了将一个或多个配体对接于单个受体结构的功能  
运行 `pyCADD` 选择模式1并使用功能6 提供包含一个多个配体的单一 `.mae` 或 `.maegz` 的文件路径即可

### Docking of One-to-Many Function

`pyCADD` 的 `Multiple Mode`提供了高性能多核并行集合式对接计算及数据处理  
请确保当前工作目录 (Current Working Directory) 在您想要保存项目文件的目录中, 并在目录中额外准备：

* 一个分行列出的, 受体蛋白所属PDB ID的列表文本文件 **.txt* 或 **.csv*
* 需要对接的外源配体 **3D结构文件**  **.pdb*  **.mae* 或其他Schrodinger支持的格式 (可选)
  
以下是一个示例 *example.txt*

    3OAP
    5JI0
    4K6I
如果晶体包含有多个不同名称的小分子配体, 请在文件中指明 (以逗号分隔) : *example.txt*

    3OAP,9CR
    5JI0,BRL
    4K6I,9RA

然后 运行 `pyCADD` 的 `Multiple Mode` , 使用功能1、2分别载入受体与配体信息, 并根据UI提示进行您需要的操作  
`pyCADD` 使用了rich库实现进度条展示, 以便您能够得知当前工作进展情况。

### 注意：

* 脚本将自动识别文本文件中每一行的PDB ID并自行批量下载  
* 应用程序默认使用系统最大核心数量进行集合式对接工作 目前尚未支持自定义核心数量  
* 所有集合式对接工作完成后, 将自动提取重要的对接结果数据, 并保存在 `result` 目录下的 `_DOCK_FINAL_RESULTS.csv`字样的文件中, 且将产生汇总文件`matrix.csv` 及 `TOTAL.csv`
* 如您仅需提取已完成工作的数据而不执行集合式对接, 请在载入受体与配体信息后使用 `Multiple Mode` 的功能8


* * *
此脚本仅限于学习和批评使用, 请勿用作其他用途。  
源码仅包含中文注释。

YH. W  
School of Pharmaceutical Sciences, Xiamen University  
2021-12-15
