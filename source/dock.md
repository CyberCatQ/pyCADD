# Dock User Guide

Dock 模块是通过串联分子对接软件 python API 来实现各类分子设计与评估功能、筛选等工作流的自动化模块。   
目前, Dock 仅针对 Schrodinger Suites 进行开发，并计划在未来添加更多主流分子对接软件的支持。

**关于Schrodinger虚拟环境**  
由于Schrodinger Suite使用了闭源环境提供Python API调用, 使用 `Dock` 模块相关功能(包括CLI命令)脚本时, 应该保证您正处于由 Schrodinger 建立的虚拟环境中，并在该环境中安装 `pyCADD`。 
* 现在，Dock 模块的 CLI 接口 `pycadd-dock` 已能够自动创建 schrodinger 虚拟 python 环境，并在其中运行 Dock 模块，您无需进行额外的`source`等操作。 `pycadd-dock` 通过 `SCHRODINGER` 变量确定您的 Schrodinger 安装位置，如果您的 `SCHRODINGER` 变量未指定，则在首次运行`pycadd-dock`接口时，会提示您提供 Schrodinger Suite 的安装目录路径。
* 通过 python 导入 Dock 模块仍然需要您自行创建或激活 Schrodinger 虚拟环境。例如使用 
    ```
    $SCHRODINGER/run schrodinger_virtualenv.py schrodinger.ve
    source schrodinger.ve/bin/activate
    python3 -m pip install pycadd -U
    ```
    命令可在当前目录下创建名为 `schrodinger.ve` 的虚拟环境，激活并在其中安装 `pyCADD`。 更多信息请参阅 [如何建立Schrodinger Python Virtual Env](https://content.schrodinger.com/Docs/r2022-1/python_api/intro.html#per-user-virtual-environments-for-installing-additional-modules)。

## Download PDB file(s)
`pyCADD` 的 `pycadd-dock download` 提供了从PDB数据库进行多线程批量下载PDB结构文件的功能（需要互联网）。
使用 `pycadd-dock download -i PDBID` 下载单一结构，或在一个文件 `example.csv` 中列出您想要下载的PDB ID：

```
1FBY
3OAP
4K6I
5JI0
```
然后使用 `pycadd-dock download -f example.csv [-s SAVE_DIR]` (在当前目录或指定目录SAVE_DIR) 批量下载它们。

## Docking of Many-to-many (Ensemble Docking)

`pyCADD` 的 `pycadd-dock ensemble-dock` 提供了高性能多核并行系综对接 (Ensemble Docking) 计算及数据处理的CLI接口。  

请确保当前工作目录在您想要保存项目文件的目录中, 并在目录中额外准备：

* 一个分行列出的, 包含所有受体蛋白所属 PDB ID的输入文件，支持的格式包括 **.csv** **.ini** **.in** **.yaml** **.yml**, 支持的输入文件格式请参阅 [为系综对接构建Dock输入文件的示例](#为系综对接构建dock输入文件的示例)。该输入文件可由`Demand`模块自动生成，或由您手动创建;  
* 一个需要对接的化合物库 **3D结构文件** **.pdb*  **.mae* **.sdf* 或其他Schrodinger支持的格式  
    如果晶体包含有多个不同名称(ID)的小分子配体, 请在文件中指明 (以逗号分隔) 以定义对接格点的中心;   
    如包含相同名称的多个小分子，则 `Dock` 会自动选择其中之一定义为对接格点的中心。

以下是一个输入文件的示例 *example.csv*：  
```
3OAP,9CR
5JI0,BRL
4K6I,9RA
```

还需要提供一个准备就绪的用于对接的化合物库文件(其中含有若干个结构)，如 *library_file.mae*，其中包含若干个已经完成了氢原子添加、质子化状态计算、能量最小化等预处理步骤的三维化合物结构。

然后 您可以通过命令
```bash
pycadd-dock ensemble-dock example.csv library_file.mae -n 24 -p SP [--del_water] [-redock] [-O]
```
启动 Ensemble Docking。  
- -n / --parallel <NUM> 指定使用的并行核心数量，如未设定将默认使用系统最大核心数量的 75%。
- -p / --precision <SP | XP | HTVS> 指定对接精度，默认为 SP。
- --del_water 指定是否删除晶体中的所有水分子。默认为保留口袋中心5埃范围内的水分子，其余删除。
- --redock 指定是否进行回顾性对接。即首先将PDB晶体的共晶小分子对接回口袋中，默认为不进行。
- -O / --overwrite 指定是否覆盖已有结果。当前路径下如已运行过ensemble docking并留有结果文件时，默认将跳过已执行过的准备工作及对接，指定该参数则重新运行完整工作流并覆盖已有文件。

Dock 模块将会下载每一个受体PDB结构，完成晶体准备过程，保留输入文件中指定小分子所结合的单一蛋白链，生成对接网格文件，就如同常规分子对接流程的蛋白准备的工作流一样。  
然后，将*library_file.mae* 中的每一个化合物交叉对接到每一个受体中。对接结束后，**M** 个受体与 **N** 个配体的 Ensemble Docking 将会得到 **M x N** 的分数矩阵。单次对接不成功时, 矩阵中的对应位置将会留空。  

使用 
```bash
pycadd-dock ensemble-dock --help
```
以获取更多帮助信息。

所有集合式对接工作完成后, 将自动提取重要的对接结果数据, 并保存在 `result` 目录下的 `_DOCK_FINAL_RESULTS.csv`字样的文件中, 且将产生汇总矩阵文件`matrix.csv`，该文件可直接用于 `Dance` 模块中用于机器学习建模等过程。

此外, `result/docking_failed_SP.csv` 文件中将会记录此轮中未能成功对接的配体-受体对。重新运行 `pycadd-dock ensemble-dock` 时，该文件中的记录配对将被自动跳过以节约时间成本；如果您希望重新对接这些配体-受体对，请删除或重命名该文件。

如果 Ensemble Docking 环节由于某些原因或意外被迫中断，您可以重新运行命令，已完成的部分将会自动跳过。您也可以使用`-O/--overwrite`参数来完全从头运行整个Ensemble Docking。  

数据提取环节意外中断时，可以使用数据提取命令来重新尝试提取及合并数据。数据提取仅尝试提取已有的结果文件，不会重新运行对接。未能成功提取的配体-受体对将会在矩阵中留空。
```bash
pycadd-dock extract-data example.csv ligand_file.mae -n 24 -p SP [--redock]
```

##  Generate Docking Reports 

CLI命令 `pycadd-dock generate-report` 提供了将若干重要分子对接于多个配体并生成对比报告的功能。  
与Ensemble Docking相似，晶体被下载并预处理完成后，所有晶体的共晶配体分子会被首先对接回晶体自身，以完成回顾性对接(Self-Docking)，并得到共晶分子的参考分数。接下来，化合物库中的分子继续依次被对接到蛋白晶体中。    
最后，将自动为化合物库中的每一个配体生成一个独立的EXCEL报告文件，以便于直观比较当前分子与共晶分子的得分情况；该分子及共晶分子的2D结构图也会被自动保存在报告文件中。  

使用命令
```bash
pycadd-dock generate-report example.csv ligand_file.mae -n 24 -p SP [--del_water] [-redock] [-O]
```
等待共晶分子的回顾性对接，及化合物库对接完成后，将会在当前目录下的 `result/` 文件夹生成对接结果文件及EXCEL报告文件。 

本质上，该接口运行与`ensemble-dock`相同的工作流，并在最后增加了EXCEL报告的生成。因此，如您已经完成了`ensemble-dock`的工作流，预计对接流程将被跳过，您可以直接使用 `generate-report` 接口来快速生成EXCEL报告。  

报告文件生成应该聚焦于少数重点分子，不建议为过多分子生成报告文件。


## 为系综对接构建Dock输入文件的示例

`pycadd-dock ensemble-dock` 支持三种输入文件格式：
* csv
* ini | in
* yaml | yml

**Note**
* 通过 `csv` 构建Dock输入文件，可以满足 `pycadd-dock ensemble-dock` 的需要，但由于其未含有受体通用名称 (如RXRα), 在调用 `pycadd-dock generate-report` 时，可能会报错，请改用另外两种格式。
* 通过 `ini` 构建Dock输入文件时，由于ini文件类型的特性，同一PDB ID不能出现在同一受体下中两次或以上。因此, 如果晶体含有多个共结晶配体，且需要分别计算，请将所有共结晶配体名称赋值于同一PDB ID下，并以英文逗号分隔。
* 通过 `yml/yaml` 构建Dock输入文件是推荐格式，可以满足所有情况下的需求，但需要注意：
    * 由于yaml语法要求，英文 `:` 后需要有一空格
    * yaml语法要求同一Section中的内容 应该具有相同的缩进，类似于python
    * yaml可以为同一晶体指定多个共结晶配体, 通过多行分割，并在每一行以 `-` 分割
    * 当共结晶配体名称仅由数字组成时，使用引号将其包裹(e.g '056')
* 可以使用 `Demand` 模块来为单个蛋白快速调研、生成 `Dock` 输入文件，仅需提供蛋白唯一识别号UniProt ID即可。请参阅 `Demand` 模块的使用说明。

以下是三种格式的示例：  

*csv*
```
1XJ7,DHT
1XQ3,R18
2AM9,TES
2AM9,DTT
2YLP,TES
2YLP,056
```

*ini | in*
```
[P10275]
    1XJ7: DHT
    1XQ3: R18
    2AM9: TES,DTT
    2YLP: TES,056
```

*yaml | yml*
```
P10275:
    1XJ7: DHT
    1XQ3: 
    - R18
    2AM9: 
    - TES
    - DTT
    2YLP:
    - TES
    - '056'
```