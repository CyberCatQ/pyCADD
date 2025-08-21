# Demand User Guide

pyCADD-Demand 是快速检索PDB数据库的工具，它可以根据用户输入的蛋白Uniprot ID，快速检索该蛋白在PDB数据库中符合条件的结构，并对检索结果进行下载、合并。  
Demand 模块需要互联网连接，但并非使用爬虫技术，而是直接通过 PDB 数据库所提供的 RESTful API 接口进行检索，同一时间下所得数据应该与您在PDB web 所获得的数据一致。因此，您无需担心数据质量。  
Demand 模块也能够为 Dock 模块开展 Ensemble Docking 快速生成必要的输入文件。

## 检索蛋白 Uniprot ID

为了收集蛋白信息，您需要在 [Uniprot](https://www.uniprot.org/) 数据库中确定所需研究的蛋白 Uniprot ID，它是蛋白的唯一识别号。  
如果有多个蛋白需要检索信息，请为其分别运行一次 Demand 模块。

## 检索 PDB 数据

### 通过 CLI 命令快速检索 PDB 数据

一旦您确定了uniprot ID, 可以通过下面的命令来快速检索当前蛋白的 PDB 数据。  

```bash
pycadd-demand [options] [UNIPROT_ID]
```
- UNIPROT_ID: 您需要检索的蛋白的 Uniprot ID。
- -f, --pdb_list_file: 查询指定csv文件中包含的PDB ID信息，而不是通过Uniprot ID查询。csv文件应该包含一个PDBID列。
- -p, --pdb_column: 指定csv文件中包含PDB ID信息的列名，默认为`PDBID`。

`Demand` 最初是为了生成 `Dock` 模块中用于 Ensemble Docking 功能输入文件而创造的，因此，下面的参数仅会影响 `Dock` 模块的输入文件生成，而不会影响 `Demand` 模块的运行。**无论如何设定可选参数，与当前蛋白相关的所有PDB数据都会被检索并保存。**

- -g, --generate: 生成Dock输入文件。默认情况下，不生成Dock输入文件。如果您需要生成Dock输入文件，必须提供Uniprot ID，并使用该参数。
- -o, --output_format: 设定生成输入文件的格式。Dock 支持 csv | in | ini | yml | yaml格式，但为了便于人类阅读，您可以选择所需的格式来输出，默认为yml。
- -c, --cutoff FLOAT: 在生成输入文件时的结构分辨率截断值，单位为埃(Angstrom)，高于截断值的结构将被过滤。默认情况下，不进行分辨率过滤。如果您需要检索更高分辨率（即具有一个更低的Cutoff）的结构，请使用该参数，并提供一个具体的浮点数。
- -m / --not_del_mutations: 在生成输入文件时不删除突变体。默认情况下，Demand 模块会删除突变的蛋白结构，只保留原始蛋白(Wild Type)的结构信息，突变体的分辨是根据PDB数据库返回信息确定的。如果您需要保留突变体的结构信息，请使用该参数。
- -e / --not_del_ignore: 在生成输入文件时不删除小分子中的**忽略分子**。**忽略分子**指通常不属于常规小分子化合物或配体的分子、原子或离子，主要为各类溶剂分子(如乙醇、DMSO等)、金属离子(镁、钾、磷酸根离子等)。您可以在 pyCADD.Demand.config 中修改这些分子的列表，或者使用该参数来保留这些分子。目前的版本中，它们包括：
    ```
    # 非配体的小分子
    IGNORE_LIG = ['EDO', 'DMS', 'IPA', 'TBY', 'ARS', 'EU', 'MG', 'IOD', 'ACT', 'CA', 'CAC', 'K', 'FMT', 'BU3', 'PGO', 'PE4', 'PO4', 'BR', 'NO3', 'BCT', 'ZN', 'SO4', 'CL', 'NA', 'AU', 'GOL', 'NI', 'YT3', 'PEG', 'PGE']
    ```


例如，如果您需要检索人类RXRα蛋白(Uniprot ID: P19793) 的PDB数据，可以使用下面的命令：

```bash
pycadd-demand P19793
```
随后将产生一个`query_data`目录，其中包含了所有检索到的PDB数据，整合完毕的PDB信息被保存于 `query_data/pdb/P19793.csv` 中，使用其他工具便可在本地进行进一步的分析和研究。  
使用`-g`参数将生成`P19793.yml`文件，仅用于作为 `pycadd-dock ensemble-dock` 的输入。

或者，您也可以自行收集需要特定查询的PDB ID，并通过一个包含了PDBID列的csv文件进行查询。例如，如果您有一个csv文件`pdb_list.csv`，其中包含了PDB ID信息，您可以使用下面的命令：

```bash
# 使用 -p 参数可指定PDBID的列名称 默认为'PDBID'
pycadd-demand -f pdbid_list.csv [-p PDBID]
```
随后将生成与Uniprot ID查询方式相同的结果目录和文件。

### 通过 python 包检索 PDB 数据

您也可以通过 python 包来检索 PDB 数据，这样可以更加灵活地使用 Demand 模块。

```python
# 导入 Demand 模块
from pyCADD import Demand

# 确定您的蛋白 Uniprot ID
uniprot_id = 'P19793'

# 创建查询器示例
query = Demand.QueryClient(uniprot_id)

# 与 Uniprot 和 PDB 数据库API交互并查询, 查询数据保存至 query/ 下
client.query()

# 清洗查询结果（仅针对Dock输入文件生成过程）:
# 去除Apo晶体   
# 去除配体未结合于目标链的晶体
# 去除非WideType晶体(optional)
# 去除非配体的小分子(e.g. DMS, optional)
# 去除分辨率高于Cutoff的晶体(optional)
client.clean_pdb_data(del_mutation=False, del_ignore=False, cutoff=None)

# 生成Dock输入文件
client.save('P19793.yml')

```