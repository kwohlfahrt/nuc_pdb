#!/usr/bin/env python3
from setuptools import setup

setup(
    name="Nuc PDB",
    version="0.0.1",
    description="A tool to generate PDB files from a nuc structure",
    packages=['nuc_pdb', 'nuc_pdb.pdb'],
    entry_points={'console_scripts': ['nuc_pdb=nuc_pdb.nuc_pdb:main']},
)
