from setuptools import setup, find_packages, Extension
import os

current_dir = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(current_dir, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

exec(open("cytocad/version.py").read())

setup(
    name='cytocad',
    version=__version__,
    packages=find_packages(),
    package_data={'cytocad.data': ['*.bed', '*.p']},
    include_package_data=True,
    ext_modules=[Extension('cytocad.bam_coverage', ['cytocad/bam_coverage.pyx'])],
    scripts=['cytocad/cytocad'],
    url='https://github.com/cytham/cytocad',
    download_url='https://github.com/cytham/cytocad/releases',
    license='gpl-3.0',
    author='Tham Cheng Yong',
    author_email='cytham@nus.edu.sg',
    description='Large copy-number variation detector with low-depth whole-genome sequencing data',
    keywords=['cytocad', 'copy number variation', 'CNV', 'whole genome sequencing', 'low depth', 'change point detection'],
    long_description=long_description,
    long_description_content_type="text/markdown",
    install_requires=['pandas>=1.1.5', 'numpy>=1.17.3', 'scipy>=1.2.1', 'pybedtools>=0.8.2', 'matplotlib>=2.2.3',
                      'ruptures>=1.1.3', 'pysam>=0.15.3', 'tagore>=1.1.0'],
    python_requires='>=3.6',
    classifiers=[
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ]
)
