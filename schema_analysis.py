#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============
Generate Figures for Protein Rarity Schema 

Written by Brian Miller
===============
"""

############################################################
# Import modules and functions

import matplotlib.pyplot as plt
import pandas as pd
import math
from schema_functions import (gauss, fitness, apply_styles, column_adjust, 
                           remove_dups, add_rel_asa, add_len, mut_dists, 
                           destab_analysis, relationship_analysis)
import xlrd         
import xlsxwriter   

# Cutoffs in kcal/mol for parsing ddG changes into stabilizing, neutral, and destabilizing 
CUTOFF_1 = .5
CUTOFF_2 = 1.0
# Minimal number of recorded mutations for a protein to be used in analysis
SAMPLE = 30

############################################################
# Average mutation distribution function from Tawfik

x_calc = []
freq1 = []
freq2 = []
freqt = []
P1 = .49
for pos in range(81):
    x = -6 + pos*.2
    x_calc.append(x)
    surf_dist = gauss(x, .54, .98, P1)
    core_dist = gauss(x, 2.05, 1.91,1-P1)
    tawfik_dist = surf_dist + core_dist
    freq1.append(surf_dist)
    freq2.append(core_dist)
    freqt.append(tawfik_dist)


############################################################
# Plot average mutation distribution for Tawfik data and mutational distribution for FireProtDB data

plot1 = plt.figure(1)
fig, ax = plt.subplots(nrows=1, ncols=1)
fig.set_size_inches(4, 3)

# Plot Tawfik average distribution with histagram
ax.plot(x_calc, freqt, color='black', linestyle='-', label='Tawfik data')

# Mutation Distribution for FireProtDB data 
# Newest version of FireProtDB can be download from https://loschmidt.chemi.muni.cz/fireprotdb
# Read data from FireProtDB, add length, and add average ddG values for duplicate mutations
col_list = ["pdb_id", "protein_name", "position", "wild_type", "mutation", "ddG", "asa", "sequence"]
mut_data_dups = pd.read_csv('Data Files/fireprotdb_results.csv', usecols=col_list)
mut_data = mut_data_dups.copy()
mut_data.dropna(subset = ["ddG"], inplace=True)  #remove rows without ddG values
remove_dups(mut_data)  #merge duplicate mutations by averaging ddG values
add_len(mut_data)
col_list.insert(5, 'len')
col_list.remove('sequence')
mut_data = mut_data[col_list]

# Write all ddG data to "mut_out [SAMPLE].xlsx"
file_name = "Data Files/mut_out " + str(SAMPLE) + ".xlsx"
writer = pd.ExcelWriter(file_name, engine = 'xlsxwriter')
mut_data.to_excel(writer, sheet_name='ddG', index=False)
column_adjust(mut_data, writer, 'ddG', 30)

# Plot ddG distribution for all mutations with color coding stabilizing, neutral, and destabilizing mutations 
ddg = mut_data.ddG      
Lower_COB = 30 - math.ceil(5*CUTOFF_1)
Upper_COB = 30 + math.ceil(5*CUTOFF_1)
n, bin, patches = ax.hist(ddg, bins=80, range=(-6,10), rwidth=0.8, 
                           density=True, align='mid', facecolor='#E9E9E9')

for index in range(Upper_COB,80):  
    plt.setp(patches[index], 'facecolor', '#7D7D7D')
for index in range(0,Lower_COB):  
    plt.setp(patches[index], 'facecolor', 'k')
cutoff_str = "{:.1f}".format(CUTOFF_1)
stab_str = "stabilizing  " + "(< -" + cutoff_str + ")"
dest_str = "destabilizing  " + "(> " + cutoff_str + ")"
plt.setp(patches[20], 'facecolor', 'k', label=stab_str)
plt.setp(patches[30], 'facecolor', '#E9E9E9', label='neutral')
plt.setp(patches[79], 'facecolor', '#7D7D7D', label=dest_str)
apply_styles(ax, -6, 10, 0, .5, 2, 1, '$\Delta\Delta$G (kcal/mol)', 'Fraction (%)', True)
ax.legend(bbox_to_anchor=(.6, .9))

# Save figure of FireProtDB ddG distribution with dpi = 500 and dpi = 150
file_name = "Figures/ddG Dist " + str(CUTOFF_1) + ".tif"
plt.savefig(file_name)



############################################################
# Plot beta-lactamase distrubtion shading regions for stabilizing, neutral, and destabilizing mutations
plot2 = plt.figure(2)
fig, ax = plt.subplots(nrows=1, ncols=1)
fig.set_size_inches(4, 3)
apply_styles(ax, -6, 10, 0, .5, 2, 1, '$\Delta\Delta$G (kcal/mol)', 'Fraction (%)', True)

blac_func = lambda x:gauss(x, .44, .96, .32) + gauss(x, 1.69, 1.77, 1-.32)

# Stabilizing region
x_calc = []
blac_dist = []
for index in range(0, Lower_COB*2+1):  
    x = -6 + index*.1
    x_calc.append(x)
    blac_dist.append(blac_func(x))
ax.fill_between(x_calc, 0, blac_dist, facecolor = 'k', label=stab_str)

# Neutral region
x_calc = []
blac_dist = []
for index in range(Lower_COB*2, Upper_COB*2+1):  
    x = -6 + index*.1
    x_calc.append(x)
    blac_dist.append(blac_func(x))
ax.fill_between(x_calc, 0, blac_dist, facecolor = '#E9E9E9', label='neutral')

# Destabilizing region
x_calc = []
blac_dist = []
for index in range(Upper_COB*2,161):  
    x = -6 + index*.1
    x_calc.append(x)
    blac_dist.append(blac_func(x))
ax.fill_between(x_calc, 0, blac_dist, facecolor = '#7D7D7D', label=dest_str)

ax.legend()

# Save figure of beta-lactamase ddG distribution with dpi = 500 and dpi = 150
file_name = "Figures/beta-lac " + str(CUTOFF_1) + ".tif"
plt.savefig(file_name)



############################################################
# Calculate and print data on destabilizing %, ASA, and length    

# Add destabilizing %, average ddG, and searchability and then save as Excel sheet
destab_analysis(mut_data, SAMPLE, CUTOFF_1, CUTOFF_2, writer)

# Remove data without asa value and add relative asa values
mut_data.dropna(subset = ["asa"], inplace=True)
add_rel_asa(mut_data)

# Write ASA data to "mut_out [SAMPLE].xlsx"
mut_data.to_excel(writer, sheet_name='ASA', index=False)
worksheet = writer.sheets['ASA']
workbook = writer.book
format1 = workbook.add_format({'num_format': '0%'})
worksheet.set_column('I:I', None, format1)
column_adjust(mut_data, writer, 'ASA', 30)

# Analyze, print, and save relationships between ddG, ASA, and length
relationship_analysis(mut_data, SAMPLE, CUTOFF_1, CUTOFF_2, writer)


############################################################
# Plot probability distibution for beta-lactamase 

plot3= plt.figure(3)

x_calc = []
blac = []
HisA = []
for pos in range(301):
    x = pos*.1
    x_calc.append(x)
    blac.append(fitness(x,.144,.036))
    HisA.append(fitness(x,.165,.065))
    
fig, ax = plt.subplots(nrows=1, ncols=1)
fig.set_size_inches(4, 3.3)

# Plot figure
apply_styles(ax, 0, 20, 0, 1, 5, 1, 'Nonsynonymous Mutations', 'Functional %', True)
plt.plot(x_calc, blac, color='black', linestyle='-')


# Save figure of beta-lactamase probability distribution
file_name = "Figures/Beta Prob Dist.tif"
plt.savefig(file_name)


############################################################
# Generate individual protein distributions for Tawfik proteins

tawfik_data = pd.read_excel('Data Files/mut_taw.xls')
mut_dists(tawfik_data, CUTOFF_1, CUTOFF_2)

# Write Tawfik data with destabilizing % to "mut_out [SAMPLE].xlsx"
tawfik_data.to_excel(writer, sheet_name='Tawfik Data', index=False)
worksheet = writer.sheets['Tawfik Data']
worksheet.set_column('I:N', None, format1)

# Adjust column widths
column_adjust(tawfik_data, writer, 'Tawfik Data', 30)

writer.save()
writer.close()

