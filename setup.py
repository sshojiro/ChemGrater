from setuptools import setup, find_packages
from os import path

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='ChemGrater',
    version='0.1',
    description='Module to fragment a molecule based on RDKit.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/sshojiro/ChemGrater',
    author='sshojiro',
    author_email='shojiro.shibayama@gmail.com',
    classifiers=[
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],
    keywords='chemoinformatics RDKit fragment',
    packages=find_packages(exclude=['data','examples','src']),
    python_requires='>=3.6',
)
