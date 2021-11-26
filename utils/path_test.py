import os

root_path = os.path.abspath(os.path.dirname(__file__)).split('pyCADD')[0]  # 项目总路径
print(root_path)
print(root_path.split('/')[-2])