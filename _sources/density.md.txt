# Density User Guide
Density 模块为自动化执行 Gaussian 相关计算的功能模块。  
该模块暂未提供CLI支持，因此，您需要直接运行 `pycadd-density` 并通过UI访问该模块。

Density 模块通过一个输入结构文件开始工作。

使用命令行  
```bash
pycadd-gauss [input_file_path]
```
来直接载入初始结构文件，并根据提示输入选项。

## Function
Density 能够根据输入结构文件的扩展名识别结构类型，然后，根据用户选择的功能，自动编写 Gaussian 输入文件，并启动计算任务。  
开始计算前，可以通过Density设定计算所用的核心与内存资源。

Density 目前支持以下功能：
- 结构优化
- 单点能计算
- 激发态吸收能与发射能计算
- 分子轨道cube文件提取