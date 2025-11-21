#!/usr/bin/env python3
"""
pyCADD简化分布式文档构建脚本
- Dock模块：使用虚拟环境 /home/yh/pyCADD/pycadd-dock.ve
- 其他模块：使用conda环境 /home/yh/miniconda3/envs/pyCADD
"""

import os
import sys
import subprocess
import shutil
import argparse
from pathlib import Path
from pyCADD import __version__

# 项目根目录
PROJECT_ROOT = Path(__file__).parent.parent
DOCS_DIR = PROJECT_ROOT / "docs"
SOURCE_DIR = DOCS_DIR / "source"
BUILD_DIR = DOCS_DIR / "build"

# 环境配置
DOCK_VENV = "/home/yh/pyCADD/pycadd-dock.ve"
CONDA_ENV = "/home/yh/miniconda3/envs/pyCADD"

# 子模块配置
MODULES = {
    "Dock": {
        "env_path": DOCK_VENV,
        "env_type": "venv",
        "rst_files": ["pyCADD.Dock.rst", "pyCADD.Dock.schrodinger.rst"],
        "md_files": ["dock.md"]
    },
    "Dynamic": {
        "env_path": CONDA_ENV,
        "env_type": "conda",
        "rst_files": ["pyCADD.Dynamic.rst"],
        "md_files": ["dynamic.md"]
    },
    "Dance": {
        "env_path": CONDA_ENV,
        "env_type": "conda",
        "rst_files": ["pyCADD.Dance.rst", "pyCADD.Dance.algorithm.rst"],
        "md_files": ["dance.md"]
    },
    "Demand": {
        "env_path": CONDA_ENV,
        "env_type": "conda",
        "rst_files": ["pyCADD.Demand.rst"],
        "md_files": ["demand.md"]
    },
    "Density": {
        "env_path": CONDA_ENV,
        "env_type": "conda",
        "rst_files": ["pyCADD.Density.rst"],
        "md_files": ["density.md"]
    }
}

def check_environments():
    """检查环境是否存在"""
    print("Checking environments...")
    
    dock_env = Path(DOCK_VENV)
    conda_env = Path(CONDA_ENV)
    
    if not dock_env.exists():
        print(f"❌ Dock virtual environment not found: {DOCK_VENV}")
        return False
    else:
        print(f"✓ Dock virtual environment found: {DOCK_VENV}")
    
    if not conda_env.exists():
        print(f"❌ Conda environment not found: {CONDA_ENV}")
        return False
    else:
        print(f"✓ Conda environment found: {CONDA_ENV}")
    
    return True

def get_sphinx_command(env_path, env_type):
    """获取对应环境的sphinx-build命令"""
    env_path = Path(env_path)
    
    if env_type == "venv":
        sphinx_build = env_path / "bin" / "sphinx-build"
    elif env_type == "conda":
        sphinx_build = env_path / "bin" / "sphinx-build"
    else:
        return "sphinx-build"
    
    if sphinx_build.exists():
        return str(sphinx_build)
    else:
        print(f"Warning: sphinx-build not found in {env_path}, using system sphinx-build")
        return "sphinx-build"

def create_module_conf(module_name):
    """为模块创建conf.py"""
    return f'''
import os
import sys
sys.path.insert(0, '../../../pyCADD')

project = 'pyCADD.{module_name}'
copyright = '2025, Yuhang Wu'
author = 'Yuhang Wu'

try:
    from pyCADD import __version__
    release = __version__
except ImportError:
    release = '0.0.1'

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.doctest',
    'sphinx.ext.intersphinx',
    'sphinx.ext.todo',
    'sphinx.ext.coverage',
    'sphinx.ext.mathjax',
    'sphinx.ext.ifconfig',
    'sphinx.ext.viewcode',
    'sphinx.ext.githubpages',
    'recommonmark',
    'sphinx_markdown_tables'
]

templates_path = ['_templates']
exclude_patterns = []
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
source_suffix = ['.rst', '.md']
master_doc = 'index'

# 自动文档配置
autodoc_default_options = {{
    'members': True,
    'member-order': 'bysource',
    'special-members': '__init__',
    'undoc-members': True,
    'exclude-members': '__weakref__'
}}
'''

def create_module_index(module_name, module_config):
    """为模块创建index.rst"""
    content = f'''
{module_name} Module Documentation
{'=' * (len(module_name) + 22)}

.. toctree::
   :maxdepth: 4
   :caption: Contents

'''
    
    # 添加文件
    for rst_file in module_config['rst_files']:
        rst_name = rst_file.replace('.rst', '')
        content += f"   {rst_name}\n"
    
    for md_file in module_config['md_files']:
        md_name = md_file.replace('.md', '')
        content += f"   {md_name}\n"
    
    content += '''

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
'''
    return content

def build_module(module_name, module_config):
    """构建单个模块的文档"""
    print(f"\n🔨 Building {module_name} module documentation...")
    print(f"   Environment: {module_config['env_type']} ({module_config['env_path']})")
    
    # 创建临时构建目录
    temp_dir = BUILD_DIR / f"{module_name}_temp"
    source_dir = temp_dir / "source"
    
    if temp_dir.exists():
        shutil.rmtree(temp_dir)
    source_dir.mkdir(parents=True, exist_ok=True)
    
    # 复制源文件
    for rst_file in module_config['rst_files']:
        src = SOURCE_DIR / rst_file
        if src.exists():
            shutil.copy2(src, source_dir)
            print(f"   Copied: {rst_file}")
    
    for md_file in module_config['md_files']:
        src = SOURCE_DIR / md_file
        if src.exists():
            shutil.copy2(src, source_dir)
            print(f"   Copied: {md_file}")
    
    # 创建配置文件
    with open(source_dir / "conf.py", "w") as f:
        f.write(create_module_conf(module_name))
    
    with open(source_dir / "index.rst", "w") as f:
        f.write(create_module_index(module_name, module_config))
    
    # 构建文档
    sphinx_cmd = get_sphinx_command(module_config['env_path'], module_config['env_type'])
    html_dir = temp_dir / "html"
    
    cmd = [
        sphinx_cmd,
        "-b", "html",
        "-E",  # 重建所有文件
        str(source_dir),
        str(html_dir)
    ]
    
    print(f"   Command: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True, cwd=str(temp_dir))
        print(f"✓ {module_name} documentation built successfully")
        if result.stdout:
            print(f"   Output: {result.stdout.strip()}")
        return html_dir
    except subprocess.CalledProcessError as e:
        print(f"❌ Error building {module_name} documentation:")
        print(f"   Error: {e.stderr}")
        if e.stdout:
            print(f"   Output: {e.stdout}")
        return None

def create_main_index():
    """创建主页面索引"""
    html_content = '''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>pyCADD Documentation</title>
    <style>
        body { 
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif; 
            margin: 0; padding: 40px; background-color: #f8f9fa; 
        }
        .container { max-width: 1200px; margin: 0 auto; background: white; padding: 40px; border-radius: 8px; box-shadow: 0 2px 10px rgba(0,0,0,0.1); }
        h1 { color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px; }
        .modules { display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 20px; margin-top: 30px; }
        .module { 
            border: 1px solid #e1e8ed; border-radius: 8px; padding: 20px; 
            transition: transform 0.2s, box-shadow 0.2s; background: white;
        }
        .module:hover { transform: translateY(-2px); box-shadow: 0 4px 20px rgba(0,0,0,0.1); }
        .module h2 { color: #34495e; margin-top: 0; }
        .module a { 
            display: inline-block; background: #3498db; color: white; 
            padding: 10px 20px; text-decoration: none; border-radius: 5px; 
            transition: background 0.2s;
        }
        .module a:hover { background: #2980b9; }
        .footer { margin-top: 40px; text-align: center; color: #7f8c8d; }
    </style>
</head>
<body>
    <div class="container">
        <h1>🧬 pyCADD Documentation</h1>
        <p>Welcome to the pyCADD (Python Computer-Aided Drug Design) documentation.</p>
        
        <div class="modules">
'''
    
    module_descriptions = {
        "Dock": "Molecular docking tools",
        "Dynamic": "Molecular dynamics simulation & analysis",
        "Dance": "Drug discovery algorithms and machine learning",
        "Demand": "Data fetcher from PDB database",
        "Density": "Quantum chemistry calculations"
    }
    
    for module_name in MODULES.keys():
        description = module_descriptions.get(module_name, f"{module_name} module")
        html_content += f'''
            <div class="module">
                <h2>{module_name}</h2>
                <p>{description}</p>
                <a href="{module_name.lower()}/{module_name.lower()}.html">User Guide</a>
                <a href="{module_name.lower()}/index.html">API Documentation</a>
            </div>
'''
    
    html_content += f'''
        </div>
        
        <div class="footer">
            <p>Generated on November 19, 2025 | pyCADD v{__version__}</p>
        </div>
    </div>
</body>
</html>'''
    
    return html_content

def merge_documentation():
    """合并所有模块文档"""
    print("\n📚 Merging module documentation...")
    
    final_dir = BUILD_DIR / "html"
    if final_dir.exists():
        shutil.rmtree(final_dir)
    final_dir.mkdir(parents=True, exist_ok=True)
    
    # 复制模块文档
    for module_name in MODULES.keys():
        temp_html = BUILD_DIR / f"{module_name}_temp" / "html"
        if temp_html.exists():
            target_dir = final_dir / module_name.lower()
            shutil.copytree(temp_html, target_dir)
            print(f"   ✓ Copied {module_name} documentation")
            # shutil.rmtree(BUILD_DIR / f"{module_name}_temp")
        else:
            print(f"   ❌ {module_name} documentation not found")
    
    # 创建主页面
    with open(final_dir / "index.html", "w") as f:
        f.write(create_main_index())
    
    print(f"✓ Documentation merged in: {final_dir}")
    return final_dir

def main():
    parser = argparse.ArgumentParser(description="Build pyCADD documentation")
    parser.add_argument("--module", choices=list(MODULES.keys()), help="Build specific module only")
    parser.add_argument("--check", action="store_true", help="Check environments only")
    args = parser.parse_args()
    
    print("🚀 pyCADD Documentation Builder")
    print("=" * 50)
    
    # 检查环境
    if not check_environments():
        print("\n❌ Environment check failed!")
        return 1
    
    if args.check:
        print("\n✓ All environments are ready!")
        return 0
    
    # 确保构建目录存在
    BUILD_DIR.mkdir(exist_ok=True)
    
    # 构建模块
    if args.module:
        print(f"\n🎯 Building single module: {args.module}")
        module_config = MODULES[args.module]
        result = build_module(args.module, module_config)
        if result:
            print(f"\n✓ {args.module} documentation ready at: {result}")
        else:
            print(f"\n❌ Failed to build {args.module} documentation")
            return 1
    else:
        print(f"\n🏗️  Building all modules...")
        built_modules = []
        
        for module_name, module_config in MODULES.items():
            result = build_module(module_name, module_config)
            if result:
                built_modules.append(module_name)
        
        if built_modules:
            final_dir = merge_documentation()
            print(f"\n🎉 Documentation build complete!")
            print(f"📖 Open: {final_dir}/index.html")
            print(f"📊 Built modules: {', '.join(built_modules)}")
        else:
            print("\n❌ No modules were built successfully!")
            return 1
    
    return 0

if __name__ == "__main__":
    sys.exit(main())