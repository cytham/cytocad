"""
A wrapper for tagore.

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


import os
# import pickle
# import cytocad
import logging
from subprocess import Popen, PIPE, STDOUT
# from cytocad.tagore import draw


# Wrapper to visualize CNVs using tagore
def tagore_wrapper(tagout, sample_name, wk_dir, ref_build, oformat):
    # Write tagore BED to file
    out_path = os.path.join(wk_dir, sample_name + '.tagore.bed')
    outwrite = open(out_path, 'w')
    _ = outwrite.write('\n'.join(tagout))
    outwrite.close()
    prefix = os.path.join(wk_dir, sample_name + '.ideo')
    process = Popen(['tagore', '-i', out_path, '-p', prefix, '-b', ref_build, '-ofmt', oformat, '-f'],
                    universal_newlines=True, stdout=PIPE, stderr=STDOUT)
    with process.stdout:
        log_subprocess(process.stdout)
    exitcode = process.wait()
    if exitcode != 0:
        logging.critical("Error: tagore failed")
        raise Exception("Error: tagore failed, see log")


def log_subprocess(out):
    for line in iter(out.readline, ''):
        if line != '\n':
            logging.debug(line.strip())


# Wrapper to visualize CNVs using tagore standalone
# def tagore_standalone(tagout, sample_name, wk_dir, ref_build):
#     # Write tagore BED to file
#     out_path = os.path.join(wk_dir, sample_name + '.tagore.bed')
#     outwrite = open(out_path, 'w')
#     _ = outwrite.write('\n'.join(tagout))
#     outwrite.close()
#     prefix = os.path.join(wk_dir, sample_name + '.ideo')
#     tag_args = TagoreArgs(out_path, prefix, ref_build)
#     data_dir = os.path.join(os.path.dirname(cytocad.__file__), 'data')
#     base_path = os.path.join(data_dir, 'base.svg.p')
#     base = open(base_path, 'rb')
#     __head__, __tail__ = pickle.load(base)
#     draw(tag_args, __head__, __tail__)


# Format tagore arguments
# class TagoreArgs:
#
#     def __init__(self, bed_path, sample_prefix, ref_build):
#         self.input = bed_path
#         self.prefix = sample_prefix
#         self.build = ref_build
