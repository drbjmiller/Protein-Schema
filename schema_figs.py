#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============
Generate Additional Figures for Protein Rarity Schema

Written by Brian Miller
===============
"""

############################################################
# Import modules and functions

import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import matplotlib.ticker as mtick
#import pandas as pd
#import math
from schema_functions import (fitness, apply_styles)
from math import (comb, log10, sqrt)
import numpy as np
from scipy.optimize import curve_fit
import scipy.integrate


############################################################
# Set parameters for protein

PROTEIN = "Beta"

if PROTEIN == "Beta":
    ALPHA = .095
    BETA = .019
    PLEN = 263
    TARGET = 150
    YAX_DIST = .12
    XMIN = 80
    XMAX = 100
    XAX_FUNSEQ = 150
    SEQMAX = 120
    ZOOM_MAX = 120
elif PROTEIN == "HisA":
    ALPHA = .0
    BETA = .09
    PLEN = 245
    TARGET = 140
    YAX_DIST = .20   
    XMIN = 30
    XMAX = 40 
    XAX_FUNSEQ = 80
    SEQMAX = 80
    ZOOM_MAX = 60
elif PROTEIN == "GFP":
    ALPHA = -.05
    BETA = .055
    PLEN = 238
    TARGET = 136
    YAX_DIST = .30   
    XMIN = 30
    XMAX = 50 
    XAX_FUNSEQ = 70
    SEQMAX = 70
    ZOOM_MAX = 70


############################################################
# Plot number of functional sequences x mutations away from optimal sequence and percent tolerated mutations

# Function for number of functional sequences with x mutations
def seq_numf(x, prot_len, func, *args):
    sequences = comb(prot_len, x) * 19**x
    if sequences > 1.7e308:
        print("Largest n: {:d}".format(x))
        return 1.7e308
    return sequences * func(x, *args)

def decline(x, a, b):
    return np.exp(-a*x-b*x**2)

blac_func = lambda x:fitness(x, ALPHA, BETA)
seqn_func = lambda x:seq_numf(x, PLEN, blac_func)
decline_func = lambda x:decline(x, ALPHA, BETA)

x_vals = range(0, XAX_FUNSEQ + 1)
seq_nums = list(map(seqn_func, x_vals))
seq_dist = [float(i)/sum(seq_nums) for i in seq_nums]
perc_tol = list(map(decline_func, x_vals))

print("Tail: {:.2e}".format(sum(seq_dist[61:140])))

fig, ax = plt.subplots(nrows=2, ncols=1)
fig.set_size_inches(5, 7)
plt.subplots_adjust(left=0.2,
                    bottom=0.07, 
                    right=0.9, 
                    top=0.95, 
                    wspace=0.4, 
                    hspace=0.2)

apply_styles(ax[0], 0, XAX_FUNSEQ, 0, YAX_DIST, 20, 10, 'n', 'Fraction (%)', False)
ax[0].yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1, decimals=0, 
                                                    symbol='%', is_latex=False))
ax[0].yaxis.set_major_locator(MultipleLocator(.02))
ax[0].yaxis.set_minor_locator(MultipleLocator(.01))
ax[0].set_title('Functional Sequences', y=.9)   
ax[0].plot(x_vals, seq_dist, color='black', linestyle='-')


apply_styles(ax[1], 0, XAX_FUNSEQ, 0, 1.0, 20, 10, 'Nonsynonymous Mutations', 'Fraction (%)', False)
ax[1].yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1, decimals=0, 
                                                    symbol='%', is_latex=False))
ax[1].yaxis.set_major_locator(MultipleLocator(.2))
ax[1].yaxis.set_minor_locator(MultipleLocator(.1))
ax[1].set_title('Tolerated Mutations', y=.9)   
ax[1].plot(x_vals, perc_tol, color='black', linestyle='-')

plt.savefig("Figures/"+ PROTEIN + " Dist.png", dpi=300)


############################################################
# Plot FSH Maps that display distribution of funcitonal sequences as heat map.


def seqdist_scaled(r, seq_dist):
    if r > SEQMAX:
        return(-1)
    elif seq_dist[r] < .001:
        return(-1)
    else: 
        return (-.95 + 8*seq_dist[r])

def plot_radial(max_rad, radii, rlabels, filename, figsize, seq_dist):
    rad = np.linspace(0, max_rad, num=max_rad+1, dtype=int)
    a = np.linspace(0, 2 * np.pi, 200)
    seq_distf = lambda x: seqdist_scaled(x, seq_dist)
    z_vals = list(map(seq_distf, rad))
    z = np.array(z_vals, ndmin=1)
    r, th = np.meshgrid(rad, a)
    z, th = np.meshgrid(z_vals, a)
    
    fig, ax = plt.subplots(nrows=1, ncols=1)
    fig.set_size_inches(figsize, figsize)
    plt.subplot(projection="polar")
     
    plt.pcolormesh(th, r, z, cmap = 'Blues', shading='gouraud')
     
    plt.grid()
    plt.thetagrids([])
    plt.rgrids(radii, labels=rlabels, fontsize=10)
    plt.savefig(filename, dpi=300)

radii = [XMIN, XMAX]
rlabels = ['', '']
plot_radial(ZOOM_MAX, radii, rlabels, "PLOSFig1.tif", 6, seq_dist)
radii = [XMIN, XMAX, SEQMAX]
rlabels = ["", str(XMAX), str(SEQMAX)]
plot_radial(TARGET, radii, rlabels, "Figures/" + PROTEIN + " " + str(TARGET) + ".png", 5.5, seq_dist)


############################################################
# Plot functional sequences on linear-log graph

r_vals = list(np.linspace(0, XAX_FUNSEQ, num=XAX_FUNSEQ+1, dtype=int))
seq_num = list(map(lambda x: log10(comb(PLEN, x)) + x*log10(19), r_vals))
perc_func = list(map(lambda x: log10(blac_func(x)), r_vals))
func_seqs = list(map(lambda x: log10(blac_func(x) * comb(PLEN,x)) + x*log10(19), r_vals))

fig, ax = plt.subplots(nrows=1, ncols=1)
fig.set_size_inches(4, 3.3)
apply_styles(ax, 0, XAX_FUNSEQ+1, -250, 320, 20, 10, 'n', 'Exponent', False)
plt.plot(r_vals, seq_num, color='black', linestyle='dotted', label="sequences")
plt.plot(r_vals, perc_func, color='green', linestyle='dashed', label="percent functional")
plt.plot(r_vals, func_seqs, color='blue', linestyle='solid', label="functual seqs")
ax.yaxis.set_major_locator(MultipleLocator(100))
ax.yaxis.set_minor_locator(MultipleLocator(50))
ax.legend()
ax.xaxis.set_label_coords(.5, .50)
plt.savefig("Figures/" + PROTEIN + " Seqs", dpi=600)


############################################################
# Curve fit Tawfik beta-lactamase data

mut_valsl = [1.6, 3.6, 5.5]
exp_dropl = [.61, .55, .39]
dxl = [1.2, 1.7, 2.2]
mut_valsh = [1.3, 2.2, 4.6]
exp_droph = [.48, .33, .32]
dxh = [.9, 1.5, 2.2]
#mut_valsm = list(np.linspace(0, 8, num=9))
#model_drop = list(map(lambda x: decline(x, .104, .019), mut_valsm))

fig, ax = plt.subplots(nrows=1, ncols=1)
fig.set_size_inches(4, 3)
apply_styles(ax, 0, 8, 0, 1, 1, .5, 'Nonsynonymous Mutations', 'Percentage Tolerated', True)
ax.yaxis.set_minor_locator(MultipleLocator(.1))

#plt.plot(mut_valsm, model_drop, color='black', linestyle='-', label="model (low)")
plt.plot(mut_valsl, exp_dropl, color='blue', linestyle='-', marker=".", label="low selection")
plt.errorbar(mut_valsl, exp_dropl, xerr=dxl, fmt='.k', elinewidth=1, capsize=2)

plt.plot(mut_valsh, exp_droph, color='green', linestyle='-', marker="*", label="high selection")
plt.errorbar(mut_valsh, exp_droph, xerr=dxh, fmt='.k', elinewidth=1, capsize=2)

plt.legend(loc=1)   # bbox_to_anchor=(1.00, 1.05))
plt.savefig("Figures/B-lac Perc", dpi=1000)


############################################################
# Curve fit flourescent protein data
if PROTEIN == "GFP":
    mut_num = list(np.linspace(2, 9, num=8))
    functional = [0.8761, 0.7302, 0.5213, 0.3159, 0.1699, 0.0801, 0.0413, 0.0154]  # Cutoff = 1000
    parameters, covariance = curve_fit(decline, mut_num, functional) #, p0=[.1, .03])
    alpha = parameters[0]
    beta = parameters[1]
    siga = sqrt(covariance[0][0])
    sigb = sqrt(covariance[1][1])
    param_str = "a = {:.3f} +/- {:.3f}   b = {:.3f} +/- {:.3f}".format(alpha, siga, beta, sigb)
    print(param_str)
    
    # Plot flourescent fitneess
    fig, ax = plt.subplots(nrows=1, ncols=1)
    fig.set_size_inches(4, 3)
    apply_styles(ax, 0, 10.2, 0, 1, 1, 1, 'n', 'Functional %', True)
    ax.yaxis.set_minor_locator(MultipleLocator(.1))
    
    mut_num.insert(0,1)
    functional.insert(0, 0.9122)
    model_drop = list(map(lambda x: decline(x, ALPHA, BETA), mut_num))
    plt.plot(mut_num, model_drop, color='black', linestyle='dotted', linewidth=1, label="model")
    plt.plot(mut_num, functional, color='blue', linestyle='-', linewidth=1, marker=".", label="experiment")
    plt.legend(loc=1)
    plt.savefig("Figures/Green Perc", dpi=1000)
