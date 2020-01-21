#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 19 18:38:46 2020

@author: usingh
"""

import setuptools
import os
import sys


#exit if python 2
if sys.version_info.major != 3:
    raise EnvironmentError("""This version of PlncPRO requires python 3.5 or higher. Please upgrade your python.\n To use PlncPRO with py2 please download older version from ccbb.jnu.ac.in/plncpro/index.html""")
    
#read description
with open("README.md", "r") as fh:
    long_description = fh.read()

#read version info
#cwd =os.path.abspath(os.path.dirname("__file__"))
#version = {}
#with open(os.path.join(cwd, "plncpro", "version.py")) as fp:
#    exec(fp.read(), version)
#version = version["__version__"]

#if version is None:
#    print("Error: version is missing. Exiting...", file=sys.stderr)
#    sys.exit(1)



setuptools.setup(
    name="plncpro",
    #version=version,
    version=open("plncpro/_version.py").readlines()[-1].split()[-1].strip("\"'"),
    author="Urminder Singh",
    author_email="usingh@iastate.edu",
    description="PlncPRO (Plant Long Non-Coding rna Prediction by Random fOrests) is a program to classify coding (mRNAs) and long non-coding transcripts (lncRNAs).",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/urmi-21/PLncPRO",
    packages=setuptools.find_packages(),
    include_package_data=True,
    package_data={#add framefinder files
            "": ["*.model", "framefinder", "plncpro_format_ff.sh"]
            },
    scripts=['plncpro/scripts/plncpro_format_ff.sh'],
    entry_points={
            'console_scripts': [
                    'plncpro = plncpro.__main__:main'
                    
                    ]
            },
    install_requires=[line.rstrip() for line in open("requirements.txt", "rt")],
    tests_require=["pytest"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: Unix",
    ],
    python_requires='>=3.5',
)