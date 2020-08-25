#!/usr/bin/env python

"""Install LoRALS"""

import os
import sys
import types
import subprocess

from datetime import datetime

if sys.version_info.major >= 3:
    import builtins
    maketrans = str.maketrans
else:
    import __builtin__ as builtins
    from string import maketrans


builtins.__LORALS_SETUP__ = True # type: bool

#   Get stuff from setuptools
from setuptools import setup
from setuptools import find_packages
from setuptools.command.install import install

#   Some basic information
NAME = "LoRALS" # type: str
AUTHOR = "Dafni Glinos" # type: str
AUTHOR_EMAIL = "dglinos@nygenome.org" # type: str
LICENSE = "" # type: str
DESCRIPTION = "" # type: str
PROJECT_URLS = { # type: Dict[str, str]
    'Documentation': '',
    'Source': '',
    'Tracker': ''
}

#   Set version information
VERSION = ""

#   Classifiers
CLASSIFIERS = [ # type: List[str, ...]
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
    'Programming Language :: Python :: 2',
    'Programming Language :: Python :: 2.7',
    'Programming Language :: Python :: 3',
    # 'Programming Language :: Python :: 3.5',
    #   Operating systems we support
    # 'Operating System :: Microsoft :: Windows',
    'Operating System :: MacOS :: MacOS X',
    'Operating System :: POSIX :: Linux',
    'Operating System :: Unix',
]

#   Specify Python version
#   We support Python 2.7 and 3.5 or higher
PYTHON_REQUIRES='>=2.7, !=3.0.*, !=3.1.*, !=3.2.*, !=3.3.*, !=3.4.*, <4' # type: str

#   Platforms
PLATFORMS = ['Linux', 'Mac OS-X', 'UNIX'] # type: List[str, ...]

#   Dependencies
INSTALL_REQUIRES = [ # type: List[str, ...]
    "pandas",
    "pybedtools",
    "pysam",
    "scipy",
]

#   Command-line scripts included in this module
import lorals.scripts as scripts
SCRIPTS = (script for script in dir(scripts) if isinstance(eval('scripts.' + script), types.FunctionType)) # type: Generator[str, ...]
SCRIPTS = (script for script in SCRIPTS if eval('scripts.%s.__module__' % script) == 'lorals.scripts') # type: Generator[str, ...]
SCRIPTS = (script for script in SCRIPTS if not script.startswith('_')) # type: Generator[str, ...]
SCRIPTS = ['%(script)s = %(pkg)s.scripts:%(script)s' % {'script': script, 'pkg': NAME.lower()} for script in SCRIPTS] # type: List[str, ...]

ENTRY_POINTS = { # type: Dict[str, List[str, ...]]
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
    entry_points=ENTRY_POINTS
)
