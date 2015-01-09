# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 11:12:23 2013

@author: Oleguer Sagarra <osagarra@ub.edu> 
"""

from setuptools import setup, find_packages
readme = open('README.txt').read()
setup(name='beepath_science',
      version='0.75',
      author='Oleguer Sagarra',
      author_email='osagarra@ub.edu',
      #license='',
      description='beepath.net scientiffic package',
      long_description=readme,
      packages=find_packages())
