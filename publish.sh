#!/bin/bash

# 设置错误时退出
set -e

# 保存当前工作目录
cwd=$(pwd)

# 切换到脚本所在目录
cd "$(dirname "$0")"

# 构建文档
echo "构建 HTML 文档..."
if ! make html; then
    echo "错误: 文档构建失败"
    exit 1
fi

# 设置路径（使用绝对路径）
docs_dir="$(cd "$(dirname "$0")" && pwd)"
build_dir="$docs_dir/html"
push_dir=$HOME/html

# 检查构建目录是否存在且不为空
if [ ! -d "$build_dir" ] || [ -z "$(ls -A "$build_dir" 2>/dev/null)" ]; then
    echo "错误: 构建目录 $build_dir 不存在或为空"
    exit 1
fi

echo "准备部署目录..."
# 如果部署目录已存在，先清空
if [ -d "$push_dir" ]; then
    echo "清空已存在的部署目录..."
    rm -rf "$push_dir"
fi

mkdir -p "$push_dir"
cd "$push_dir"

# 初始化 Git 仓库
echo "初始化 Git 仓库..."
git init

# 设置用户信息（如果没有设置的话）
if ! git config user.name >/dev/null 2>&1; then
    git config user.name "GitHub Actions"
fi
if ! git config user.email >/dev/null 2>&1; then
    git config user.email "actions@github.com"
fi

# 检查并设置远程仓库
echo "设置远程仓库..."
if git remote get-url origin >/dev/null 2>&1; then
    # 远程仓库已存在，更新 URL
    git remote set-url origin https://github.com/CyberCatQ/pyCADD.git
else
    # 远程仓库不存在，添加新的
    git remote add origin https://github.com/CyberCatQ/pyCADD.git
fi

# 获取远程分支信息
echo "获取远程分支信息..."
git fetch origin

# 检查并切换到 gh-pages 分支
echo "切换到 gh-pages 分支..."
if git ls-remote --heads origin gh-pages | grep -q gh-pages; then
    # 远程 gh-pages 分支存在
    echo "远程 gh-pages 分支已存在，切换到该分支..."
    git checkout -b gh-pages origin/gh-pages
    # 清空当前分支的所有文件
    git rm -rf . 2>/dev/null || true
else
    # 远程 gh-pages 分支不存在，创建新分支
    echo "创建新的 gh-pages 分支..."
    git checkout --orphan gh-pages
fi

# 复制构建的文件
echo "复制构建文件..."
cp -r "$build_dir"/* .

# 添加 .nojekyll 文件（GitHub Pages 优化）
touch .nojekyll

# 提交更改
echo "提交更改..."
git add .
if git diff --staged --quiet; then
    echo "没有更改需要提交"
else
    git commit -m "Deploy to GitHub Pages - $(date '+%Y-%m-%d %H:%M:%S')"
    
    # 推送到远程
    echo "推送到 GitHub Pages..."
    git push origin gh-pages
    echo "部署完成！"
fi

rm -rf "$push_dir"
# 返回原目录
cd "$cwd"