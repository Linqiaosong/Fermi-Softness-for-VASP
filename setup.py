# coding: utf-8

from setuptools import setup, find_packages

setup(
    name='FermiSoftness',
    version='1.2.0.1',
    author='linqiaosong',
    author_email='linqiaosong@outlook.com',
    url='https://github.com/Linqiaosong/Fermi-Softness-for-VASP',
    description=u'A script for calculating Fermi-Softness.',
    packages=find_packages(),
    install_requires=[
        'numpy>=1.22.0',
        'ase>=3.22.1'       
        ],
    entry_points={
        'console_scripts': [ ]},
    python_requires='>=3.8'
    )