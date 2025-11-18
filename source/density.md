# Density User Guide
Density 模块为自动化执行 Gaussian 相关计算的功能模块。  
2.0版本开始，该模块已提供CLI支持，因此，您可直接运行 `pycadd-density` 命令来提交所需的分子结构计算任务。

使用命令行  
```bash
pycadd-density --help
pycadd-density optimize --help
pycadd-density single-point --help
pycadd-density resp --help
pycadd-density resp2 --help
```
查看可用命令及选项。


Density 目前支持以下功能：
- 结构优化
- 单点能计算
- RESP/RESP2 电荷计算


## 参数说明

### 结构优化
```bash
pycadd-density optimize [-c charge] [-m multiplicity] [-d dft] [-b basis_set] [-s solvent] [-l] [-p cpu_num] [-M memory] [-S save_dir] STRUCTURE_FILE
```
- `-c charge`：分子电荷，默认为0。
- `-m multiplicity`：分子多重度，默认为1。
- `-d dft`：所使用的DFT泛函，默认为B3LYP。
- `-b basis_set`：所使用的基组，默认为6-31G(d)。
- `-s solvent`：所使用的溶剂模型，默认为无溶剂（真空）。使用时指定溶剂名称，如water、methanol等。
- `-l, --loose`：使用宽松的收敛标准进行优化。
- `-p cpu_num`：使用的CPU核心数，默认为1。
- `-M memory`：分配给Gaussian计算的内存，默认为4GB。
- `-S save_dir`：计算结果保存目录，默认为当前目录。
- `STRUCTURE_FILE`：输入的分子结构文件，支持格式包括XYZ、SDF等。

### 单点能计算
```bash
pycadd-density single-point [-c charge] [-m multiplicity] [-d dft] [-b basis_set] [-s solvent] [-e] [-p cpu_num] [-M memory] [-S save_dir] STRUCTURE_FILE
```
- `-c charge`：分子电荷，默认为0。
- `-m multiplicity`：分子多重度，默认为1。
- `-d dft`：所使用的DFT泛函，默认为B3LYP。
- `-b basis_set`：所使用的基组，默认为6-31G(d)。
- `-s solvent`：所使用的溶剂模型，默认为无溶剂（真空）。使用时指定溶剂名称，如water、methanol等。
- `-p cpu_num`：使用的CPU核心数，默认为1。
- `-e, --esp_calculate`：计算分子静电势。
- `-M memory`：分配给Gaussian计算的内存，默认为4GB。
- `-S save_dir`：计算结果保存目录，默认为当前目录。
- `STRUCTURE_FILE`：输入的分子结构文件，支持格式包括XYZ、SDF等。

### RESP/RESP2 电荷计算
```bash
pycadd-density resp [-c charge] [-m multiplicity] [-d dft] [-b basis_set] [-s solvent] [-p cpu_num] [-M memory] [-S save_dir] STRUCTURE_FILE

pycadd-density resp2 [-c charge] [-m multiplicity] [-d dft] [-b basis_set] [-s solvent] [-D delta] [-p cpu_num] [-M memory] [-S save_dir] STRUCTURE_FILE
```
- `-c charge`：分子电荷，默认为0。
- `-m multiplicity`：分子多重度，默认为1。
- `-d dft`：所使用的DFT泛函，默认为B3LYP。
- `-b basis_set`：所使用的基组，默认为6-31G(d)。
- `-s solvent`：所使用的溶剂模型，默认为无溶剂（真空）。使用时指定溶剂名称，如water、methanol等。
- `-D delta`：RESP2方法中的delta参数，默认为0.5，仅适用于resp2命令。
- `-p cpu_num`：使用的CPU核心数，默认为1。
- `-M memory`：分配给Gaussian计算的内存，默认为4GB。
- `-S save_dir`：计算结果保存目录，默认为当前目录。
- `STRUCTURE_FILE`：输入的分子结构文件，支持格式包括XYZ、SDF等。