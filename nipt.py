import pandas as pd
import numpy as np
import os
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import sys
import vcf
import os
import re
import json
import math
import six
from matplotlib import colors
from itertools import cycle
from colormath.color_objects import LCHabColor, sRGBColor
from colormath.color_conversions import convert_color
from scipy.stats import gaussian_kde
from sklearn import svm

coverage = pd.read_table('cov_female.bins', names=['chrom', 'start', 'end', 'count', 'nzpos', 'length', 'nzratio'])

def apply_dropped_spine(ax, spines=('left', 'bottom'), xgrid=False, smart_bounds=False,
                        drop=5):
    if drop:
        for sp in spines:
            ax.spines[sp].set_position(('outward', drop))

    for loc, spine in ax.spines.items():
        if loc in spines:
            spine.set_smart_bounds(smart_bounds)
            spine.set_linewidth(0.8)
        else:
            spine.set_color('none') # don't draw spine

    # turn off ticks where there is no spine
    if 'left' in spines:
        ax.yaxis.set_ticks_position('left')
    elif 'right' in spines:
        ax.yaxis.set_ticks_position('right')
    else:
        # no yaxis ticks
        ax.yaxis.set_ticks([])

    if 'bottom' in spines:
        ax.xaxis.set_ticks_position('bottom')
    elif 'top' in spines:
        ax.xaxis.set_ticks_position('top')
    else:
        # no xaxis ticks
        ax.xaxis.set_ticks([])

    ax.yaxis.grid(True, linestyle='-', alpha=0.1, linewidth=1.5)
    if xgrid:
        ax.xaxis.grid(True, linestyle='-', alpha=0.1, linewidth=1.5)
    else:
        ax.xaxis.grid(False)

    ax.tick_params('both', which='major', width=0.8, direction='out', pad=3)

    leg = ax.get_legend()
    if leg is not None:
        ltext  = leg.get_texts()
        llines = leg.get_lines()
        frame  = leg.get_frame()

        #frame.set_facecolor('0.80')
        frame.set_linewidth(0.8)
        plt.setp(ltext, fontsize='12')
        plt.setp(llines, linewidth=1.5)

chrom_human_order = lambda x: 'chr{:03d}'.format(int(x[3:])) if x[3:].isdigit() else x

binsize = coverage['length'].value_counts().index[0]

covs = coverage[(coverage['length'] == binsize) & (coverage['chrom'].apply(len) <= 6)].copy()
covs['chrom_sortkey'] = covs['chrom'].apply(chrom_human_order)
covs = covs.sort_values(by=['chrom_sortkey', 'start']).reset_index(drop=True)
print(covs.head())

def colormap_lch(n, lum=(35, 65), chroma=75, start=0, end=300):
    if isinstance(lum, list) or isinstance(lum, tuple):
        lum = cycle(lum)
    else:
        lum = cycle([lum])

    rgbs = []
    for hue, lumn in zip(np.linspace(start, end, n), lum):
        rgb = convert_color(LCHabColor(lumn, chroma, hue), sRGBColor).get_value_tuple()
        rgbs.append('#{:02x}{:02x}{:02x}'.format(*map(int, np.array(rgb).clip(0, 1) * 255.0)))
    return rgbs

plt.style.use('ggplot')
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = 'Ubuntu'
plt.rcParams['font.monospace'] = 'Ubuntu Mono'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 10
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.titlesize'] = 10
plt.rcParams['xtick.labelsize'] = 8
plt.rcParams['ytick.labelsize'] = 8
plt.rcParams['legend.fontsize'] = 10
plt.rcParams['figure.titlesize'] = 12

plt.figure(figsize=(9, 3))

chromosomes = covs['chrom'].unique()
plt.xticks(np.arange(1, len(chromosomes) + 1),
           [c[3:] for c in chromosomes])

plt.xlabel('chromosome')
plt.ylabel('log2 read count/bin')

colors_ = list(six.iteritems(colors.cnames))

for name, rgb in six.iteritems(colors.ColorConverter.colors):
    hex_ = colors.rgb2hex(rgb)
    colors_.append((name, hex_))

# Transform to hex color values.
hex_ = [color[1] for color in colors_]
# Get the rgb equivalent.
rgb = [colors.hex2color(color) for color in hex_]
# Get the hsv equivalent.
hsv = [colors.rgb_to_hsv(color) for color in rgb]

chroms = sorted(covs['chrom_sortkey'].unique())
chromcolors = {chrom: color for chrom, color
               in zip(chroms, hex_)}
chromcolors['chrY'] = '#3d3d3d'

chroms = sorted(covs['chrom_sortkey'].unique())
chromcolors = {chrom: color for chrom, color
               in zip(chroms, colormap_lch(len(chroms), end=360))}
chromcolors['chrY'] = '#3d3d3d'

fig, ax = plt.subplots(1, 1, figsize=(20, 8))

chromnameY = 8.5
wiggle_y = .1

for i, (chrom, rows) in enumerate(covs.groupby('chrom_sortkey')):
    ax.scatter(rows.index, np.log2(rows['count']),
               edgecolor='none', s=4, alpha=0.2,
               c=chromcolors[chrom])

    center_y = np.median(np.log2(rows['count']).dropna())
    ax.plot([rows.index[0], rows.index[-1]],
            [center_y, center_y],
            c='black',linewidth=0.4)
    
    center_x = np.mean(rows.index)
    ax.annotate(chrom[3:].lstrip('0'),
                (center_x, chromnameY),
                color=chromcolors[chrom], ha='center')

apply_dropped_spine(ax)

ax.set_ylabel('Reads per 100k bin (log2)')
ax.set_xlim(-100, 15800*2)
ax.set_ylim(4, 16)

fig.savefig('nipt.png')

