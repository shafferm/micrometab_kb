from setuptools import setup, find_packages
import requests
from os import path

__author__ = "shafferm"


setup(
    name='micrometab_kb',
    version='0.1',
    install_requires=["requests", "flask", "py2cytoscape", "sqlalchemy", "networkx", "biom-format"],
    packages=find_packages(),
    url='https://github.com/shafferm/micrometab_KB/',
    license='BSD',
    author='Michael Shaffer',
    author_email='michael.shaffer@ucdenver.edu',
    description='A web service for metabolic exploration of bacteria based on 16S identity.',
)
