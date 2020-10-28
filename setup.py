#!/usr/bin/env python3

"""Install LoRALS"""

import os
import sys
import glob
import types
import builtins

from datetime import datetime
from typing import Dict, Iterator, List


builtins.__LORALS_SETUP__: bool = True

#   Get stuff from setuptools
from setuptools import setup
from setuptools import find_packages
from setuptools.command.install import install

#   Some basic information
NAME: str = "LoRALS"
AUTHOR: str = "Dafni Glinos"
AUTHOR_EMAIL: str = "dglinos@nygenome.org"
LICENSE: str = ""
DESCRIPTION: str = ""
PROJECT_URLS: Dict[str, str] = {
    'Documentation': '',
    'Source': '',
    'Tracker': ''
}

#   Set version information
VERSION = ""

#   Classifiers
CLASSIFIERS: List[str] = [
    #   See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    #   For more classfiiers
    #   How mature is this project?
    'Development Status :: 3 - Alpha',
    #   What environment does this run in?
    'Environment :: Console',
    #   Indicate who your project is intended for
    'Intended Audience :: Science/Research',
    'Intended Audience :: End Users/Desktop',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    # 'Topic :: Software Development :: Libraries :: Python Modules',
    #   Language
    'Natural Language :: English',
    #   Pick your license as you wish (should match "license" above)
    # 'License :: OSI Approved :: MIT License',
    #   Specify the Python versions you support here. In particular, ensure
    #   that you indicate whether you support Python 2, Python 3 or both.
    'Programming Language :: Python :: 3 :: Only',
    #   Operating systems we support
    # 'Operating System :: Microsoft :: Windows',
    'Operating System :: MacOS :: MacOS X',
    'Operating System :: POSIX :: Linux',
    'Operating System :: Unix',
    'Typing :: Typed',
]

#   Specify Python version
#   We support Python 3.5 or higher
PYTHON_REQUIRES: str ='>=3.5, <4'

#   Platforms
PLATFORMS: List[str] = ['Linux', 'Mac OS-X', 'UNIX']

#   Dependencies
INSTALL_REQUIRES: List[str] = [
    "pandas",
    "pybedtools",
    "pysam",
    "scipy",
]

#   Command-line scripts included in this module
import lorals.scripts as scripts

SCRIPTS: Iterator[str] = (script for script in dir(scripts) if isinstance(eval('scripts.' + script), types.FunctionType))
SCRIPTS: Iterator[str] = (script for script in SCRIPTS if eval('scripts.%s.__module__' % script) == 'lorals.scripts')
SCRIPTS: Iterator[str] = (script for script in SCRIPTS if not script.startswith('_'))
SCRIPTS: List[str] = ['%(script)s = %(pkg)s.scripts:%(script)s' % {'script': script, 'pkg': NAME.lower()} for script in SCRIPTS]

ENTRY_POINTS: Dict[str, List[str]] = {
    'console_scripts': SCRIPTS
}

#   Run setup
setup(
    name=NAME,
    version=VERSION,
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    description=DESCRIPTION,
    license=LICENSE,
    project_urls=PROJECT_URLS,
    platforms=PLATFORMS,
    python_requires=PYTHON_REQUIRES,
    classifiers=CLASSIFIERS,
    install_requires=INSTALL_REQUIRES,
    packages=find_packages(),
    cmdclass={'install': install},
    package_data={"lorals": ['blacklists/*.bed', 'py.typed']},
    entry_points=ENTRY_POINTS,
    scripts=glob.glob('scripts/*.sh')
)
