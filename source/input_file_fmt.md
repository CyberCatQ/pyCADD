## 构建Dock输入文件时的规范

pyCADD支持三种输入文件格式：
* csv
* ini | in
* yaml | yml


### Note
* 通过 `csv` 构建Dock输入文件，可以满足 `pycadd-dock ensemble-dock` 的需要，但由于其不含有受体通用名称 (如RXRα) , 在调用 `pycadd-dock quick-report` 时，可能会报错，请改用另外两种格式。
* 通过 `ini` 构建Dock输入文件时，由于ini文件类型的特性，同一PDB ID不能出现在同一受体下中两次或以上。因此, 如果晶体含有多个共结晶配体，且需要分别计算，请将所有共结晶配体名称赋值于同一PDB ID下，并以英文逗号分隔。
* 通过 `yaml` 构建Dock输入文件是推荐格式，可以满足所有情况下的需求，但需要注意：
    * 由于yaml语法要求，英文 `":"` 后需要有一空格
    * yaml语法要求同一Section中的内容 应该具有相同的缩进，类似于python
    * yaml可以为同一晶体指定多个共结晶配体, 通过多行分割，并在每一行以 `"-"` 分割
    * 当共结晶配体名称仅由数字组成时，使用引号将其包裹(e.g '056')

以下是三种格式的示例：  
*csv*

    1XJ7,DHT
    1XQ3,R18
    2AM9,TES
    2AM9,DTT
    2YLP,TES
    2YLP,056

*ini | in*

    [P10275]
        1XJ7: DHT
        1XQ3: R18
        2AM9: TES,DTT
        2YLP: TES,056


*yaml | yml*

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