#!/usr/bin/env python

import os
import sys

from distutils.core import setup

scripts = ['scripts/cherry_sight.py',
          ]

# setup and install
setup(name='adipo_sight',
      version='0.5',
      author='Adam Labadorf',
      author_email='alabadorf@gmail.com',
      package_dir={'':'src'},
      py_modules=['chipsequtil.db','chipsequtil.mww','chipsequtil.log',]
      packages=['adipo_sight'],
      package_data={'': ['org_settings.cfg']},
      scripts=scripts,
     )
