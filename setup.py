from setuptools import setup, find_packages
from setuptools.command.build_py import build_py
import glob
import os
import pkg_resources

from snipit import __version__, _program

setup(name='snipit',
      version=__version__,
      packages=find_packages(),
      scripts=["snipit/scripts/snp_functions.py"],
      install_requires=[
            "biopython>=1.70",
            "matplotlib>=3.2.1"
        ],
      description='snipit',
      url='https://github.com/aineniamh/snipit',
      author='Aine OToole',
      author_email='aine.otoole@ed.ac.uk',
      entry_points="""
      [console_scripts]
      {program} = snipit.command:main
      """.format(program = _program),
      include_package_data=True,
      keywords=[],
      zip_safe=False)
