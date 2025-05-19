from setuptools import setup
from os import path

__copyright__ = "Copyright 2022-2025, The PSMPA Project."


def get_version(rel_path):
    '''get version information from __init__.py file'''
    with open(rel_path) as fp:
        cont = fp.read().splitlines()
    for line in cont:
        if line.startswith('__version__'):
            delim = '"' if '"' in line else "'"
            return line.split(delim)[1]
    else:
        raise RuntimeError("Unable to find version infomation.")


curr_dir = path.dirname(path.abspath(__file__))
version = get_version(path.join(curr_dir, 'psmpa', '__init__.py'))

long_description = (
    "PSMPA is a Python pipeline to predict secondary metabolism potential"
    " using 16S/18S rRNA or amplicans for a single strain or microbial communities.")

setup(name='psmpa',
      version=version,
      license="GPL",
      description=('PSMPA: Prediction of Secondary Metabolism Potential using Amplicans'),
      author='Zhen-Yi Zhou',
      author_email="gavinchou64@gmail.com",
      url='https://github.com/BioGavin/psmpa',
      classifiers=["Programming Language :: Python :: 3.8",
                   "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
                   "Operating System :: OS Independent"],
      packages=['psmpa'],
      scripts=['scripts/psmpa1', 'scripts/psmpa2', 'scripts/psmpa-fungi'],
      include_package_data=True,
      python_requires=">=3.8",
      install_requires=['numpy', 'biopython', 'pandas', 'biom-format', 'joblib'],
      package_data={'psmpa': ['Rscripts/*.R',
                              'default_files/psmpa1/*.gz',
                              'default_files/psmpa1/pro_ref/*',
                              'default_files/psmpa2/*.tsv.gz',
                              'default_files/psmpa2/blast_db/*',
                              'default_files/psmpa_fungi/*.gz',
                              'default_files/psmpa_fungi/blast_db/*']},
      long_description=long_description)
