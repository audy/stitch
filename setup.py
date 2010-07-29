try:
    from setuptools import setup
except:
    from distutils import *
import os

DOC = 'See README.md'

setup(
    name="stitch",
    version="1.0.0",
    author="Austin G. Davis-Richardson",
    author_email="harekrishna@gmail.com",
    description="Assembler for Illumina Overlapping Paired-End Reads",
    url="http://www.github.com/audy/stitch",
    license="GPLv3",
    long_description=DOC,
    
    install_requires = [],
    
    packages = ["stitch"],
    
    entry_points = {
        'console_scripts': [ 
            'stitch = stitch.stitch:main',
            ]
    }

)
