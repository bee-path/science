# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 11:12:23 2013

@author: Oleguer Sagarra <osagarra@ub.edu> 
"""

from setuptools import setup, find_packages
readme = open('README.md').read()
setup(name='beepath_science',
      version='0.78',
      author='Oleguer Sagarra & Mario Gutierrez-Roig',
      author_email='osagarra@ub.edu',
      license='GPLv3',
      description="This python package contains the class definition and functions used in the analysis of the data from the Bee-path experiment.\n It may be used to analyze point-like GPS general mobility data and to aggregate such points into displacements and paused states.\nFor more details check out the webpage of the project  [Bee-path](http://bee-path.net/?lang=en) and the related publications.\nThe module is fully documented in a *quasi-standard*, *pythonic* way (but could be improved).",
      long_description=readme,
      packages=find_packages(),
      include_package_data = True
      )
