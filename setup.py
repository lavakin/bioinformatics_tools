#!/usr/bin/env python3
from setuptools import setup

setup(
    name='bioinf',
    version='0.0.1',
    description='Bioinformatics tools',
    packages=['bioinf'],
    install_requires=[
        'biopython',
        'fastaparser',
        'numpy',
        'matplotlib',
        'scipy'
    ],
    include_package_data=True,
    zip_safe=False,
)
