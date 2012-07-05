#!/usr/bin/env python

import os
import sys

from distutils.core import setup

scripts = ['scripts/cherry_sight.py',
          ]

# these files do not exist in the repo, must be downloaded
# and expanded separately!
data_files = ['data/adipo_sight.db',
              'data/motif_names.txt',
             ]
data_dirs = ['tmpl',
             'images',
            ]

# setup and install
setup(name='adipo_sight',
      version='0.5',
      author='Adam Labadorf',
      author_email='alabadorf@gmail.com',
      package_dir={'':'src'},
      py_modules=['adipo_sight.db','adipo_sight.mww','adipo_sight.log',],
      packages=['adipo_sight'],
      package_data={'': data_files,
                    'data': data_dirs},
      scripts=scripts,
     )
