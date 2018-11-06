#!/usr/bin/env python3
"""
Setup tools installation.

Run this script to install or upgrade multiprot.
"""

import os
from setuptools import setup, find_packages

REQUIREMENTS = open(os.path.join(os.path.dirname(__file__),
                                 'requirements.txt')).readlines()
setup(name='multiprot',
      version='1.0.0',
      install_requires=REQUIREMENTS,
      packages=find_packages(),
      scripts=['multiprot/scripts/multipr'],
      author="Francisco Javier Guzman-Vega",
      description=("Automated pipeline to build protein models connecting two \
         or more structured domains with disordered linkers."),
      author_email="francisco.guzmanvega@kaust.edu.sa",
      url="https://github.com/StruBE-KAUST/multiprot",
      python_requires='>=3',
      )
