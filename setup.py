#!/usr/bin/env python

import os
import sys

from distutils.core import setup

scripts = ['scripts/cherry_sight.py',
          ]

data_dirs = ['data/tmpl',
            ]

# setup and install
setup(name='adipo_sight',
      version='0.5',
      author='Adam Labadorf',
      author_email='alabadorf@gmail.com',
      package_dir={'adipo_sight':'src'},
      py_modules=['adipo_sight.db','adipo_sight.mww','adipo_sight.log',]
      packages=['adipo_sight'],
      package_data={'adipo_sight': data_dirs},
      scripts=scripts,
     )
