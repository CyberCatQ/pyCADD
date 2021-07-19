automatedMD
==========
尝试自动化调用Schrodinger Python API 执行晶体准备与对接，以及AMBER分子动力学模拟准备、运行和MD轨迹的基本分析

## Platform  
* Linux  
## Required    
* [Schrodinger Suite](https://www.schrodinger.com/)2018-2 或更高版本
* [AMBER](http://ambermd.org/) 18 或更高版本(要使用PMEMD 请安装CUDA 9.0或以上版本)
* [Gaussian](http://gaussian.com/) 09.D01 或更高版本 (脚本默认使用g16)
* [Multiwfn 3.7](http://sobereva.com/multiwfn/) (用以计算RESP2电荷)

### Python Version
* 3.7 or Higher


### Python Modules  
* pytraj  
* xlsxwriter  
* Openbabel

如果缺少python库 请使用`pip install <package>`命令安装

## Script Function
|        Name        | Function |
| -----------------  | -------- |
|*py4schrodinger.py* | 自动执行PDB晶体获取、优化、格点文件生成、对接等命令(请使用$SCHRODINGER/run运行此脚本) |
|*pyMDprepare.py*    | 准备AMBER MD必要文件(拓扑及坐标、力场参数等文件)   
|*RESP2.sh*          | 调用Gaussian执行坐标优化并计算RESP2(0.5)电荷(已集成于 *pyMDprepare.py*)  |
|*pyPMEMD.py*        | 调用AMBER PMEMD(GPU加速)执行能量最小化、体系加热与分子动力学模拟  |
|*pynalysis.py*      | MD轨迹分析 输出RMSD/RMSF、氢键、二面角等变化情况 提取最低势能构象  调用MMPBSA计算吉布斯自由能变/熵变  |
|*libpynalysis/*        | pynalysis.py所需模块包|

## How to Use
使用`git clone`或`cp`命令 将 `automatedMD/` 文件夹移入需要分析的PDB晶体文件夹(以其ID命名)的目录下  
(以下简称根目录)，不需要额外文件   
然后，您可以从 `bin/` 文件夹中执行需要的应用。

以PDBID 1XLS为例：  

    mkdir 1XLS
    cp -r automatedMD/ 1XLS/
    cd 1XLS/automatedMD/bin
    chmod +x ./*

    #根据需要运行
    ./py4schrodinger
    ./pyMDprepare  
    ./runPMEMD   
    ./pynalysis  

更多帮助信息，请参阅API文档。

## Schrodinger Multiple Docking
除直接执行 `bin/py4schrodinger` 并使用UI操作外， 脚本可还进行自动一对多、多对一对接工作。  

### Docking of Many-to-One Function
请将包含多个配体的单一`.mae`或`.maegz`文件与脚本置于被命名为 *PDBID* 的根文件夹内，按照一对一对接相同的操作：直接运行`py4schrodinger` 并选择 功能6 , 按照提示输入文件名即可。

### Docking of One-to-Many Function
请将此脚本文件至于一个单独的文件夹中, 并在文件夹中额外准备：

* 一个分行列出的，受体蛋白所属PDB ID的列表文本文件 *.txt* 
* 需要对接的外源配体 3D结构文件  *.pdb* & *.mae*或其他Schrodinger支持的格式 (可选)
  
以下是一个示例 *example.txt*

    3OAP
    5JI0
    4K6I
如果晶体包含有多个不同名称的小分子配体，请在文件中指明，并以逗号分割(还可特别指明构象类型): *example.txt*

    3OAP,9CR,agonist
    5JI0,BRL,agonist
    4K6I,9RA,

然后 请使用额外参数运行`py4schrodinger`：  

    ./py4schrodinger -r|--receptor <receptors list file> [-l|--ligand <ligand file> -p|--precision <precision> -n|--cpu <cpus> -k|--no-check ]

参数含义如下：

* -r|--receptor  <受体蛋白所属PDB ID的列表文本文件路径>
* -l|--ligand <外源配体文件路径(可选)>
* -p|--precision <对接精确度>
* -n|--cpu <允许同时运行的进程数/使用的CPU核心数>
* -k|--no-check 无需确认直接启动对接进程
  
当没有-l参数传入时，脚本仅执行受体列表内的晶体获取、处理与内源配体对接，而不会将任何外源配体对接到列表中的受体。  

也可通过

    ./py4schrodinger -h  

命令获取详情。  

请注意：

* 脚本将自动识别.txt文本文件中每一行的PDB ID并自行批量下载、处理蛋白并总是执行内源配体对接 计算对接结果的RMSD值  
* 此后，脚本按照指定的精度将外源配体对接到列表中的所有受体上(如果有)
* 为了高效的完成对接任务 脚本在处理含有多条链的蛋白结构时，将仅保留**单链**（同源多聚体）或指定的内源配体所在**单链**（异源多聚体）进行对接工作
* 所有对接工作完成后，脚本还将自动提取重要的对接结果数据，并保存在带有`_FINAL_RESULTS.csv`字样的文件中。


### NOTE
* 如需进行外源性小分子与靶点蛋白的分子动力学模拟 请将小分子、受体蛋白PDB格式文件依次命名为`xxxlig_YOURLIGNAME.pdb` `xxxpro_YOURLIGNAME.pdb`  (xxx为命名的任意代号 至少三位字母或数字, YOURLIGNAME为自定义的配体名称)  
并移至脚本所在文件夹根目录下 在提示输入PDB ID时手动输入xxx即可

* 为了顺利运行所有脚本功能  
请按照 `py4schrodinger` → `pyMDprepare` → `runPMEMD` → `pynalysis` 依次运行脚本。
* * *
此脚本仅限于学习和批评使用，请勿用作其他用途。   
源码仅包含中文注释。

YH. W  
School of Pharmaceutical Sciences, Xiamen University  
2021-07-13  