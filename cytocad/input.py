"""
Functions for parsing and verifying input parameters.

Copyright (C) 2021 Tham Cheng Yong

This file is part of CytoCAD.

CytoCAD is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

CytoCAD is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with CytoCAD.  If not, see <https://www.gnu.org/licenses/>.
"""


import sys
import argparse
from cytocad import __version__


# Parse input
def input_parser(args=sys.argv[1:]):
    parser = argparse.ArgumentParser(description="CytoCAD is a tool for discovering large genomic copy-number variation \
through coverage anormaly detection (CAD) using low-depth whole-genome sequencing data.",
                                     formatter_class=argparse.RawTextHelpFormatter, usage=msg())

    parser.add_argument("input", type=str,
                        metavar="[BAM]",
                        help="""path to mapped BAM file.
Format: .bam""")

    parser.add_argument("dir", type=str,
                        metavar="[work_directory]",
                        help="""path to work directory. Directory will be created 
if it does not exist.""")

    parser.add_argument("-b", "--build", type=str, metavar="str",
                        default='hg38',
                        help="""build version of human reference genome assembly [hg38]""")

    parser.add_argument("-c", "--colors", nargs='+', metavar="hex_color",
                        default=None,
                        help="""hex color for neutral, gain, and loss CNVs on chromosome
ideograms respectively separated by space ['#a6a6a6' '#990000' '#000099']""")

    parser.add_argument('-f', '--format', type=str, metavar='[png/pdf]',
                        default='png',
                        help="Output format of chromosome illustration figure [png]")

    parser.add_argument("-i", "--interval", type=int, metavar="int",
                        default=50000,
                        help="""spread between each point in a chromosome where "
coverage is enquired, in bp. Minimum CNV sensitive 
detection size ~= interval*rolling [50000]""")

    parser.add_argument("-j", "--buffer", type=int, metavar="int",
                        default=10,
                        help="buffer window size of each point, in bp [10]")

    parser.add_argument("-r", "--rolling", type=int, metavar="int",
                        default=10,
                        help="rolling mean window size [10]")

    parser.add_argument("-p", "--penalty", type=int, metavar="int",
                        default=500,
                        help="""Linear kernel penalty value for change 
point detection using Ruptures [500]""")

    parser.add_argument("-s", "--scale", type=float, metavar="float",
                        default=0.25,
                        help="""proportion of mean coverage to be used for 
buffering to call hetero- and homozygous CNVs 
(E.g. a heterozygous loss is where a coverage (c) 
satisfies: mean-mean*scale <= c < mean+mean*scale [0.25]""")

    parser.add_argument("--add_plots", action='store_true',
                        help="output additional coverage plots in 'fig' directory")

    parser.add_argument("--debug", action='store_true',
                        help="run in debug mode")

    parser.add_argument("-v", "--version", action='version',
                        version=__version__,
                        help="show version and exit")

    parser.add_argument("-q", "--quiet", action='store_true',
                        help="hide verbose")

    args = parser.parse_args(args)
    return args


# Custom usage message
def msg():
    return "cytocad [options] [BAM] [WORK_DIRECTORY]"
