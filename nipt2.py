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
from rpy2 import robjects as ro
from sklearn import svm

coverage = pd.read_table('cov.bins', names=['chrom', 'start', 'end', 'count', 'nzpos', 'length', 'nzratio'])

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

nuccomp = pd.read_table("gc.txt",usecols=range(5))
nuccomp.columns = ['chrom', 'start', 'end', 'pct_at', 'pct_gc']
print(nuccomp.head())

nuccomp['pct_gc'] = nuccomp['pct_gc'] / (nuccomp['pct_gc'] + nuccomp['pct_at'])
del nuccomp['pct_at']

print(nuccomp.head())
print(covs)
print(nuccomp)
covs_with_gc = pd.merge(covs, nuccomp, left_on=['chrom', 'start', 'end'],right_on=['chrom', 'start', 'end']).dropna()
covs_with_gc.plot(kind='scatter', x='pct_gc', y='count',edgecolor='none', s=4, c='black')

print(covs_with_gc.head())

gc_count_mtx = covs_with_gc[['pct_gc', 'count']].as_matrix().T
density = gaussian_kde(gc_count_mtx, bw_method=0.1)(gc_count_mtx)

fig, axes = plt.subplots(1, 2, figsize=(10, 4))

covs_with_gc.plot(kind='scatter', x='pct_gc', y='count', edgecolor='none', s=5,
                  c=density, vmin=0, cmap='jet', ax=axes[0])
covs_with_gc.plot(kind='scatter', x='pct_gc', y='count', edgecolor='none', s=5,
                  c=density, vmin=0, cmap='jet', ax=axes[1])
axes[1].set_ylim(-1, 20000)
axes[1].set_xlim(0.3, 0.65)

plt.tight_layout()

def loess_fit(x, y, px=None, model=None, alpha=0.5):
    if model is None:
        model = ro.r('y ~ x')

    if px is None:
        px = np.linspace(min(x), max(y), 22)[1:-1]

    fitframe = ro.DataFrame({'x': ro.FloatVector(x), 'y': ro.FloatVector(y)})
    loessmodel = ro.r.loess(model, fitframe, span=alpha)

    predframe = ro.DataFrame({'x': ro.FloatVector(px)})
    predy = ro.r.predict(loessmodel, predframe)
    preddata = [(x, predy[i]) for i, x in enumerate(px)]

    return np.array(preddata).transpose()

gc_count_mtx_training = gc_count_mtx[:, gc_count_mtx[1] > 90]
model = svm.OneClassSVM(nu=0.06, kernel="rbf", gamma=0.5).fit(gc_count_mtx_training.T)
isinlier = model.predict(gc_count_mtx.T) == 1

gc_count_inliers = gc_count_mtx[:, isinlier]
fit_x, fit_y = loess_fit(gc_count_inliers[0], gc_count_inliers[1],
                               np.linspace(0, 1, 200), alpha=.4)

fig, axes = plt.subplots(1, 2, figsize=(10, 4))

covs_with_gc.plot(kind='scatter', x='pct_gc', y='count', edgecolor='none', s=5,
                  c=isinlier, vmin=-.1, vmax=1.1, cmap='RdBu_r', ax=axes[0])
covs_with_gc.plot(kind='scatter', x='pct_gc', y='count', edgecolor='none', s=5,
                  c=isinlier, vmin=-.1, vmax=1.1, cmap='RdBu_r', ax=axes[1])
axes[1].set_ylim(-1, 20000)
axes[1].set_xlim(0.3, 0.65)

axes[1].plot(fit_x, fit_y, c='magenta', lw=3)
axes[0].plot(fit_x, fit_y, c='magenta', lw=3)

plt.tight_layout()

gc_inliers_center = np.vstack(loess_fit(gc_count_inliers[0], gc_count_inliers[1],
                                              gc_count_inliers[0], alpha=.4))

gc_logenr = gc_count_inliers.copy()
gc_logenr[1] = np.log2(gc_count_inliers[1] / gc_inliers_center[1])

orig_density = gaussian_kde(gc_count_inliers, bw_method=0.1)(gc_count_inliers)
adj_density = gaussian_kde(gc_logenr, bw_method=0.1)(gc_logenr)

fig, axes = plt.subplots(1, 2, figsize=(10, 4))

axes[0].scatter(gc_count_inliers[0], gc_count_inliers[1], c=orig_density, edgecolor='none', s=5)
axes[0].set_xlabel('pct_gc')
axes[0].set_ylabel('mapped reads')
axes[0].set_ylim(-1, 20000)
axes[0].set_xlim(0.3, 0.65)

axes[1].scatter(gc_logenr[0], gc_logenr[1], c=adj_density, edgecolor='none', s=5)
axes[1].set_xlabel('pct_gc')
axes[1].set_ylabel('adjusted mapped reads')
axes[1].set_ylim(-1, 1)

plt.tight_layout()

_, expected = loess_fit(gc_count_inliers[0], gc_count_inliers[1],
                              covs_with_gc['pct_gc'], alpha=.4)
covs_w_logenr = covs_with_gc.copy()
covs_w_logenr['log2enrichment'] = np.log2(covs_w_logenr['count'] / expected)
covs_w_logenr = covs_w_logenr[covs_w_logenr['count'] > 0].dropna()

chrY_enrichment = covs_w_logenr[covs_w_logenr['chrom'] == 'chrY']['log2enrichment']
xv = np.linspace(-8, 4, 100)
yv = gaussian_kde(chrY_enrichment, bw_method=0.2)(xv)
plt.plot(xv, yv)
plt.xlabel('log2 enrichment')
plt.ylabel('density')
plt.title('chromosome Y')
xv[np.argmax(yv)]

fig, ax = plt.subplots(1, 1, figsize=(10, 5))

chromnameY = -1.2
wiggle_y = .05

for i, (chrom, rows) in enumerate(covs_w_logenr.groupby('chrom_sortkey')):
    ax.scatter(rows.index, rows['log2enrichment'],
               edgecolor='none', s=4, alpha=0.4,
               c=chromcolors[chrom])

    center_y = np.median(rows['log2enrichment'].dropna())
    ax.plot([rows.index[0], rows.index[-1]],
            [center_y, center_y],
            c='black',linewidth=0.4)
    
    center_x = np.mean(rows.index)
    ax.annotate(chrom[3:].lstrip('0'),
                (center_x, chromnameY + (wiggle_y * (i % 2) * 2 - wiggle_y)),
                color=chromcolors[chrom], ha='center')

apply_dropped_spine(ax)

ax.set_ylabel('Reads per 100k bin (log2)')
ax.set_xlim(-100, 31000)
ax.set_ylim(-4, 1)

fig.savefig('nipt_recal.png')
plt.close(fig)
