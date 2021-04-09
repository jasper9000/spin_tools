#!/usr/bin/env python

#from distutils.core import setup
import setuptools

setuptools.setup(name='spin_tools',
      version='1.0',
      description='Spin Tools',
      author='Jasper Riebesehl',
      author_email='jariebes@physnet.uni-hamburg.de',
      package_dir={"": "spin_tools"},
      packages=setuptools.find_packages(where="spin_tools"),
     )