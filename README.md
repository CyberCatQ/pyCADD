automatedMD
==========
尝试自动化调用Schrodinger Python API 执行晶体准备与对接，以及AMBER分子动力学模拟准备、运行和MD轨迹的基本分析

## Platform  
* Linux  
## Required    
* Schrodier Suite 2018-2 或更高版本
* AMBER 18 或更高版本(要使用PMEMD 请安装CUDA 9.0或以上版本)
* Gaussian 09.D01 或更高版本 
* Multiwfn 3.7
### Python Version
* 3.5.2 or Higher
### Python Packages  
* pytraj  
* xlsxwriter  
  
## Script Function
|        Name        | Function |
| -----------------  | -------- |
|*py4schrodinger.py* | 自动执行PDB晶体获取、优化、格点文件生成、对接等命令(请使用$schrodinger/run运行此脚本) |
|*pyMDprepare.py*    | 准备AMBER MD必要文件(拓扑及坐标、力场参数等文件)   
|*RESP2.sh*          | 调用Gaussian执行坐标优化并计算RESP2(0.5)电荷(已集成于 *pyMDprepare.py*)  |
|*pyPMEMD.py*        | 调用AMBER PMEMD(GPU加速)执行能量最小化、体系加热与分子动力学模拟  |
|*runPMEMD.sh*       | 选择要运行MD的GPU并启动 pyPMEMD.py   |
|*pynalysis.py*      | MD轨迹分析 输出RMSD/RMSF、氢键、二面角等变化情况 提取最低势能构象  调用MMPBSA计算吉布斯自由能变/熵变  |
|*pynalysis/*        | pynalysis.py所需模块包|

## How to Use
将 `automatedMD/` 内所有脚本文件及`pynalysis/`文件夹移入需要分析的PDB晶体文件夹(以其ID命名)根目录下  
(不需要额外文件)

以PDBID 1XLS为例：  

    mkdir 1XLS
    cp -r automatedMD/* 1XLS/
    cd 1XLS/
    chmod +x *

    #根据需要运行
    run py4schrodinger.py   #请使用schrodinger run运行py4schrodinger.py
    ./pyMDprepare.py
    ./runPMEMD.sh   
    ./pynalysis.py
### NOTE
* 如需进行外源性小分子与靶点蛋白的分子动力学模拟 请将小分子、受体蛋白及二者复合物PDB格式文件依次命名为`xxxlig.pdb` `xxxpro.pdb` `xxxcom.pdb`  (xxx为命名的任意代号 至少三位字母或数字)  
并移至脚本所在文件夹根目录下 在提示输入PDB ID时手动输入xxx即可

* 为了顺利运行所有脚本功能  
请按照 `py4schrodinger.py` → `pyMDprepare.py` → `runPMEMD.sh` → `pynalysis.py` 依次运行脚本。
* * *
此脚本仅限于学习和批评使用，请勿用作其他用途。   
源码仅包含中文注释。
