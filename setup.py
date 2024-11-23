'''
Date: 2024-05-06 17:18:10
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2024-11-18 21:57:54
Description: 
'''

import json
import pathlib
import sys

# Always prefer setuptools over distutils
from setuptools import find_packages, setup

here = pathlib.Path(__file__).parent.resolve()

# Get the long description from the README file
long_description = (here / "README.md").read_text(encoding="utf-8")
requirements = json.loads((here / "requirements.json").read_text())
version_info = (here / "lazydock/__version__.py").read_text(encoding="utf-8")
for line in version_info.split('\n'):
    if '__version__' in line:
        __version__ = line[line.find('"')+1:-1]
    if '__author_email__' in line:
        __author_email__ = line[line.find('"')+1:-1]
    if '__author__' in line:
        __author__ = line[line.find('"')+1:-1]
    if '__url__' in line:
        __url__ = line[line.find('"')+1:-1]
        
# decide which dyna lib to compare with
if sys.platform == "win32":
    dynlib = [
    ]
elif sys.platform in ["linux", "linux2"]:
    dynlib = [
    ]
else:
    dynlib = [
    ]
    

setup(
    name = "lazydock",
    version = __version__,

    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Utilities",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        # "Programming Language :: Python :: 3 :: Only",
    ],
        
    keywords = ["molecular docking", "Utilities"],
    description = "A Python package for molecular docking",
    long_description = long_description,
    long_description_content_type='text/markdown',
    python_requires=">=3.7, <3.13",
    license = "MIT Licence",

    url = __url__,
    author = __author__,
    author_email = __author_email__,
    
    packages = find_packages(exclude=["test", "test."]),
    include_package_data = True, # define in MANIFEST.in file
    package_data = {"lazydock": dynlib},
    
    entry_points={
        "console_scripts": [
            "lazydock-pml=lazydock_pymol_plugin.main:exec_from_cli",
            "lazydock-cli=lazydock.scripts:main"
        ],
    },
    
    platforms = "any",
    
    install_requires=requirements['std'],
)

# pip install .

# python setup.py sdist
# twine upload dist/lazydock-0.6.2.tar.gz