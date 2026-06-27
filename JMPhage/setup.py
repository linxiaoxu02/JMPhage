from setuptools import setup, find_packages

setup(
    name="jmphage",
    version="1.0.0",
    author="LXX, YMS",
    description="A modular pipeline for phage analysis",
    # 核心改动 1：自动寻找项目下的所有包（包括 modules）
    packages=find_packages(), 
    
    # 核心改动 2：确保非 Python 文件（如可能存在的配置文件）也能被打包
    include_package_data=True,
    
    # 核心改动 3：入口点映射到新的包路径
    # 格式：命令名 = 包名.文件名:函数名
    entry_points={
        'console_scripts': [
            'jmphage = jmphage.jmphage:main', 
        ],
    },
    
    # 核心改动 4：指定 Python 版本要求（可选，建议加上）
    python_requires='>=3.8',
    zip_safe=False,
)