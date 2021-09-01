"""
Detect coverage changes by change point detection

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
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pybedtools import BedTool
from collections import OrderedDict, defaultdict
import ruptures as rpt
import cytocad
from cytocad.depth_limit import normalbed


# Coverage anomaly detection
def cad(
        subdata,
        mean_cov,
        upper_cov,
        sample_name,
        ref_build='hg38',
        interval=50000,
        interval_buf=10,
        rolling_size=10,
        penalty=500,  # 500
        zygo_scale=0.25,
        cov_plots=False,
        wk_dir='./',
        colors=None  # [Neutral, Gain, Loss]
):
    # Define colors
    if colors is None:
        neutral_color, gain_color, loss_color = ['#a6a6a6', '#990000', '#000099']
    else:
        if len(colors) == 3:
            neutral_color, gain_color, loss_color = colors
        else:
            raise Exception('Error: The variable "colors" has to be a list of length 3, indicating the hex color of Neutral, '
                            'Gain and Loss CNVs in this sequential order.')

    # Define data directory
    data_dir = os.path.join(os.path.dirname(cytocad.__file__), 'data')

    # Define file paths according to reference build
    if ref_build == 'hg38':
        filter_path = os.path.join(data_dir, 'hg38_curated_filter_main_arte.bed')
        ideo_path = os.path.join(data_dir, 'hg38_ucsc_ideogram.bed')
        main_chr_path = os.path.join(data_dir, 'hg38_sizes_main.bed')
    else:
        raise Exception("Error: Reference genome build %s is not recognised. CytoCAD only supports build hg38." % ref_build)

    # Define other variables
    buf = mean_cov * zygo_scale
    het_loss = mean_cov / 2
    hom_loss = 0
    het_gain = mean_cov + mean_cov / 2
    hom_gain = mean_cov * 2

    # Make filtered-region BedTool object
    filter_bed = BedTool(filter_path)

    # Create intervals of genomic coordinates for each chromosome
    main_chr = set()
    chr_range = OrderedDict()
    with open(main_chr_path) as f:
        for line in f:
            chrm, start, end = line.split('\t')
            chr_range[chrm] = range(int(start), int(end), interval)
            main_chr.add(chrm)
    sort_chr = list(main_chr)
    sort_chr.sort()

    # Make read alignment BedTool object
    bed_str = normalbed(subdata)
    bed = BedTool(bed_str, from_string=True)
    bed = bed.sort()

    # For each chromosome, analyse read coverage by BED intersect
    data = defaultdict(list)
    region = defaultdict(list)
    chrx_avg = 0
    chry_avg = 0
    for chromo in chr_range:
        # print('Processing ' + chromo)
        cov_dict = {}
        region_dict = {}
        string = ''
        for i in chr_range[chromo]:
            string += chromo + '\t' + str(max(i - interval_buf, 1)) + '\t' + str(i + interval_buf) + '\n'
            cov_dict[chromo + '-' + str(max(i - interval_buf, 1))] = 0
            region_dict[chromo + '-' + str(max(i - interval_buf, 1))] = max(i, 1)
        interval_bed = BedTool(string, from_string=True)
        interval_bed = interval_bed.sort()
        intersect1 = interval_bed.intersect(filter_bed, wa=True, v=True)  # Remove gap regions, prevent false amplifications
        intersect2 = intersect1.intersect(bed, wa=True, wb=True)  # Intersect with read alignment data
        intersect3 = intersect1.intersect(bed, wa=True, v=True)  # Not intersected with read alignment data
        # Get number of read coverage at each interval
        cov_list = []
        for line in intersect2:
            line = list(line)
            cov_list.append(line[0] + '-' + line[1])
        for key in cov_list:
            cov_dict[key] += 1
        if chromo == 'chrX':
            total = 0
            nocov = 0
            for key in set(cov_list):
                total += cov_dict[key]
            for _ in intersect3:
                nocov += 1
            chrx_avg = total / (len(set(cov_list)) + nocov)
            # print('Chrx avg: ' + str(chrx_avg))
        elif chromo == 'chrY':
            total = 0
            nocov = 0
            for key in set(cov_list):
                total += cov_dict[key]
            for _ in intersect3:
                nocov += 1
            chry_avg = total / (len(set(cov_list)) + nocov)
            # print('Chry avg: ' + str(chry_avg))
        # For each region falling in filter bed, assign mean coverage
        intersect4 = interval_bed.intersect(filter_bed, wa=True)
        for line in intersect4:
            line = list(line)
            if chromo == 'chrX':
                cov_dict[line[0] + '-' + line[1]] = chrx_avg
            elif chromo == 'chrY':
                cov_dict[line[0] + '-' + line[1]] = chry_avg
            else:
                cov_dict[line[0] + '-' + line[1]] = mean_cov
        for key in cov_dict:
            data[key.split('-')[0]].append(min(cov_dict[key], upper_cov))  # coverage counts
            region[key.split('-')[0]].append(region_dict[key]/1000000)  # genomic coordinate

    # Process ideogram coordinates
    ideo_dict = {}
    with open(ideo_path) as f:
        for i in f:
            if not i.startswith('#'):
                if i.split('\t')[0] not in ideo_dict:
                    ideo_dict[i.split('\t')[0]] = {}
                ideo_dict[i.split('\t')[0]][str(i.split('\t')[1]) + '-' + str(i.split('\t')[2])] = i.split('\t')[3]

    # Process by eight chromosomes cycle
    groups = [[0, 8], [8, 16], [16, 24]]
    tagout = ['#chr\tstart\tstop\tfeature\tsize\tcolor\tchrCopy']
    out = []
    cycle = 1
    for g in groups:
        if cov_plots:
            fig = plt.figure()
        n = 1
        for chromo in sort_chr[g[0]:g[1]]:
            d = {'x': region[chromo], 'y': data[chromo]}
            df = pd.DataFrame(d)
            xcoord = df.x.to_numpy()
            rolling_mean = df.y.rolling(window=rolling_size).mean()
            signal = np.array(rolling_mean)
            signal[np.isnan(signal)] = mean_cov
            # Scale coverage for consistent change detection
            if chromo == 'chrX':
                signal_scaled = signal/chrx_avg * 8
            elif chromo == 'chrY':
                signal_scaled = signal/chry_avg * 8
            else:
                signal_scaled = signal/mean_cov * 8
            algo = rpt.KernelCPD(kernel="linear", min_size=2).fit(signal_scaled)
            result = algo.predict(pen=penalty)
            _span = []
            # if len(result) % 2 != 0:
            _result = [0] + result
            for i in range(len(_result) - 1):
                _span.append([_result[i], _result[i + 1]])
            # else:
            #     _result = [0] + result
            #     for i in range(len(_result) - 1):
            #         _span.append([_result[i], _result[i + 1]])
            ideo = {}
            cov = {}
            coord = {}
            for i in _span:
                span_ideo = []
                left = region[chromo][i[0]] * 1000000
                right = region[chromo][i[1] - 1] * 1000000
                for j in ideo_dict[chromo]:
                    if int(j.split('-')[0]) < left <= int(j.split('-')[1]):
                        span_ideo.append(ideo_dict[chromo][j])
                    if int(j.split('-')[0]) < right <= int(j.split('-')[1]):
                        span_ideo.append(ideo_dict[chromo][j])
                        break
                ideo[str(i[0]) + '-' + str(i[1])] = span_ideo
                size = len(range(i[0], i[1] - 1))
                s = 0
                for d in range(i[0], i[1] - 1):
                    s += data[chromo][d]
                cov[str(i[0]) + '-' + str(i[1])] = round(s / size, 3)
                coord[str(i[0]) + '-' + str(i[1])] = [int(left), int(right)]
            label = {}
            if chromo == 'chrX':
                if mean_cov - buf <= chrx_avg:  # XX
                    for i in cov:
                        if het_loss - buf <= cov[i] < het_loss + buf:
                            label[i] = 'loss-hetero'
                        elif hom_loss <= cov[i] < hom_loss + buf:
                            label[i] = 'loss-homo'
                        elif het_gain - buf <= cov[i] < het_gain + buf:
                            label[i] = 'gain-hetero'
                        elif hom_gain - buf <= cov[i]:
                            label[i] = 'gain-homo'
                        else:
                            label[i] = 'neutral-double'
                else:  # X or no X
                    for i in cov:
                        if hom_loss <= cov[i] < hom_loss + buf:
                            label[i] = 'loss-single'
                        elif mean_cov / 2 + buf <= cov[i]:
                            label[i] = 'gain-single'
                        else:
                            label[i] = 'neutral-single'
            elif chromo == 'chrY':
                if mean_cov / 2 - buf <= chry_avg:  # Y
                    for i in cov:
                        if hom_loss <= cov[i] < hom_loss + buf:
                            label[i] = 'loss-single'
                        elif mean_cov / 2 + buf <= cov[i]:
                            label[i] = 'gain-single'
                        else:
                            label[i] = 'neutral-single'
                else:  # no Y
                    for i in cov:
                        label[i] = 'nil-nil'
            else:
                for i in cov:
                    if het_loss - buf <= cov[i] < het_loss + buf:
                        label[i] = 'loss-hetero'
                    elif hom_loss <= cov[i] < hom_loss + buf:
                        label[i] = 'loss-homo'
                    elif het_gain - buf <= cov[i] < het_gain + buf:
                        label[i] = 'gain-hetero'
                    elif hom_gain - buf <= cov[i]:
                        label[i] = 'gain-homo'
                    else:
                        label[i] = 'neutral-double'
            for i in cov:
                copy = label[i].split('-')[0]
                zygo = label[i].split('-')[1]
                if copy == 'neutral':
                    if zygo == 'double':
                        tagout.append(
                            chromo + '\t' + str(coord[i][0]) + '\t' + str(coord[i][1]) + '\t0\t1\t' + neutral_color + '\t1'
                        )
                        tagout.append(
                            chromo + '\t' + str(coord[i][0]) + '\t' + str(coord[i][1]) + '\t0\t1\t' + neutral_color + '\t2'
                        )
                    elif zygo == 'single':
                        tagout.append(
                            chromo + '\t' + str(coord[i][0]) + '\t' + str(coord[i][1]) + '\t0\t1\t' + neutral_color + '\t1'
                        )
                elif copy == 'nil':
                    continue
                else:
                    if copy == 'gain':
                        color = gain_color
                    else:  # loss
                        color = loss_color
                    if zygo == 'hetero':
                        tagout.append(
                            chromo + '\t' + str(coord[i][0]) + '\t' + str(coord[i][1]) + '\t0\t1\t' + color + '\t1'
                        )
                        tagout.append(
                            chromo + '\t' + str(coord[i][0]) + '\t' + str(coord[i][1]) + '\t0\t1\t' + neutral_color + '\t2'
                        )
                    elif zygo == 'homo':
                        tagout.append(
                            chromo + '\t' + str(coord[i][0]) + '\t' + str(coord[i][1]) + '\t0\t1\t' + color + '\t1'
                        )
                        tagout.append(
                            chromo + '\t' + str(coord[i][0]) + '\t' + str(coord[i][1]) + '\t0\t1\t' + color + '\t2'
                        )
                    elif zygo == 'single':
                        tagout.append(
                            chromo + '\t' + str(coord[i][0]) + '\t' + str(coord[i][1]) + '\t0\t1\t' + color + '\t1'
                        )
                    out.append(
                        chromo + '\t' + str(coord[i][0]) + '\t' + str(coord[i][1]) + '\t' + ideo[i][0] + ';' + ideo[i][1] + '\t' +
                        str(round(cov[i], 1)) + '\t' + copy + '\t' + zygo
                    )
            if cov_plots:
                span = []
                if len(result) % 2 != 0:
                    result.append(0)
                    result_iter = iter(result)
                    for i in result_iter:
                        span.append([i, next(result_iter)])
                    _ = span.pop(-1)
                else:
                    result_iter = iter(result)
                    for i in result_iter:
                        span.append([i, next(result_iter)])
                ax = fig.add_subplot(2, 4, n)
                ax.plot(xcoord, signal, label='SMA', color='red', alpha=0.8)
                for p in span:
                    ax.axvspan(region[chromo][p[0] - 1], region[chromo][p[1] - 1], alpha=0.3, color='blue')
                ax.set_title(chromo)
                ax.set_ylim(0, upper_cov + 2)
                n += 1
        if cov_plots:
            n = 1
            fig.add_subplot(111, frame_on=False)
            plt.tick_params(labelcolor="none", bottom=False, left=False)
            plt.xlabel('Coordinate (MB)')
            plt.ylabel('Coverage')
            plt.tight_layout()
            fig_out_path = os.path.join(wk_dir, 'fig', sample_name + '_cov' + str(cycle) + '.svg')
            plt.savefig(fig_out_path,
                        dpi=100
                        )
            cycle += 1
    return out, tagout
