import os


src_path = os.path.dirname(__file__)
print(src_path)

if not os.path.exists(src_path + '/RESP2.sh'):  # 检查必要文件与依赖程序
    raise RuntimeError('Script RESP2.sh is not Found in %s!' % src_path)