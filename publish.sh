#!/bin/bash

# 设置错误时退出
set -e

# 保存当前工作目录
cwd=$(pwd)

# 切换到脚本所在目录
cd "$(dirname "$0")"
rm source/*.rst -f
make api

# 构建文档
echo "构建文档"
if ! make distributed; then
    echo "错误: 文档构建失败"
    exit 1
fi

# 设置路径（使用绝对路径）
docs_dir="$(cd "$(dirname "$0")" && pwd)"
build_dir="$docs_dir/build/html"
push_dir=$docs_dir/../gh-pages

# 检查构建目录是否存在且不为空
if [ ! -d "$build_dir" ] || [ -z "$(ls -A "$build_dir" 2>/dev/null)" ]; then
    echo "错误: 构建目录 $build_dir 不存在或为空"
    exit 1
fi

echo "准备部署目录..."
# 如果部署目录已存在，先清空
if [ -d "$push_dir" ]; then
    echo "清空已存在的部署目录..."
    rm -rf "$push_dir/d* $push_dir/index.html"
fi

mkdir -p "$push_dir"
cd "$push_dir"

# 获取远程分支信息
echo "获取远程分支信息..."
git fetch origin

# 复制构建的文件
echo "复制构建文件..."
cp -r "$build_dir"/* .

# 添加 .nojekyll 文件（GitHub Pages 优化）
touch .nojekyll

# 提交更改
# echo "提交更改..."
# git add .
# if git diff --staged --quiet; then
#     echo "没有更改需要提交"
# else
#     git commit -m "Update documentation"
#     # 推送到远程
#     echo "推送到 GitHub Pages..."
#     git push origin gh-pages
#     echo "部署完成！"
# fi

# 返回原目录
cd "$cwd"