
## CytoCAD - Copy-number variation caller using low-depth whole-genome sequencing data
[![Build Status](https://travis-ci.org/cytham/cytocad.svg?branch=master)](https://travis-ci.com/cytham/cytocad)
[![PyPI pyversions](https://img.shields.io/pypi/pyversions/cytocad)](https://pypi.org/project/cytocad/)
[![PyPI versions](https://img.shields.io/pypi/v/cytocad)](https://pypi.org/project/cytocad/)
[![Conda](https://img.shields.io/conda/v/bioconda/cytocad)](https://anaconda.org/bioconda/cytocad)
[![Github release](https://img.shields.io/github/v/release/cytham/cytocad?include_prereleases)](../../releases)
[![PyPI license](https://img.shields.io/pypi/l/cytocad)](./LICENSE.txt)
  
CytoCAD is a bioinformatics tool for the identification of large genomic copy-number variation through coverage anomaly detection
 (CAD) using mapped whole-genome sequencing (WGS) data. It has been tested in low-depth (~8X) Oxford Nanopore WGS long-read
  data. Its output displays chromosome illustrations demarcating regions of copy-number gains (Red) or losses (Blue).

### Basic information:
* Takes as input a mapped whole-genome sequencing BAM file and output a BED file and a chromosome ideogram figure
* Uses [Ruptures](https://github.com/deepcharles/ruptures) python package for change point detection of read coverage data per
 chromosome 
* Uses [tagore](https://github.com/jordanlab/tagore) for chromosome ideogram illustrations

## Getting Started

### Quick run

```
cytocad [Options] sample.bam working_dir 
```

| Argument | Comment |
| :--- | :--- |
| sample.bam | Input mapped WGS BAM file |
| working_dir | Working directory |

#### Output
| Output file | Comment |
| :--- | :--- |
| ${sample}.ideo.svg | Chromosome ideogram produced by [tagore](https://github.com/jordanlab/tagore) |
| ${sample}.CNV.bed | BED file of chromosome regions with CNV |

For more information, see [wiki](https://github.com/cytham/cytocad/wiki).

### Operating system: 
* Linux (x86_64 architecture, tested in Ubuntu 16.04)  

### Installation:
There are three ways to install CytoCAD:
#### Option 1: Conda (Recommended)
```
# Installing from bioconda automatically installs all dependencies 
conda install -c bioconda cytocad
```
#### Option 2: PyPI (See dependencies below)
```
# Installing from PyPI requires own installation of dependencies, see below
pip install cytocad
```
#### Option 3: GitHub (See dependencies below)
```
# Installing from GitHub requires own installation of dependencies, see below
git clone https://github.com/cytham/cytocad.git 
cd cytocad
pip install .
```

### Installation of dependencies
* bedtools >=2.26.0
* samtools >=1.3.0
* rsvg-convert >=2.40.13

Please make sure each executable binary is in PATH.
##### 1. _bedtools_
Please visit [here](https://bedtools.readthedocs.io/en/latest/content/installation.html) for instructions to install.

##### 2. _samtools_
Please visit [here](http://www.htslib.org/download/) for instructions to install.

##### 3. _rsvg-convert_
```
sudo apt-get update
sudo apt-get install librsvg2-bin
```

## Versioning
See [CHANGELOG](./CHANGELOG.txt)

## Citation

Not available yet

## Author

* **Tham Cheng Yong** - [cytham](https://github.com/cytham)

## License

This project is licensed under GNU General Public License - see [LICENSE.txt](./LICENSE.txt) for details.

## Limitations
* Chromosome pairs illustrated by tagore may resemble sister chromatids, but they are in fact homologous pairs
* Phasing of CNVs for each chromosome homologous pair is not yet possible.
* The default minimum size of detectable CNV is about 500 kb. It can be adjusted by the 'interval' and 'rolling' parameters
 following the equation: minimum size ~= interval*rolling 
* Other chromosomal structural variations, such as inversions, have to be detected by other tools, such as [NanoVar](https://github.com/cytham/nanovar). NanoVar has incorporated CytoCAD in its pipeline from version 1.4.0 onwards.
