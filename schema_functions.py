#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============
Functions uses in figures.py

Written by Brian Miller
===============
"""

import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import matplotlib.ticker as mtick
import numpy as np
import pandas as pd
import math
import statistics as st
import scipy.integrate
import xlrd
import xlsxwriter


############################################################
# Gaussian and fitness functions

def gauss(x, mean, sig, scale):
    CG = scale/math.sqrt(6.28*sig**2)
    gauss = np.exp(-(x-mean)**2/(2*sig**2))
    return CG*gauss

def fitness(x, a, b):
    fitness = np.exp(-a*x-b*x**2)
    return fitness


############################################################
# Apply styles to plots

def apply_styles(sp_axis, x1, y1, x2, y2, tick_l, tick_s, x_lab, y_lab, percent):
    ax = sp_axis
    plt.style.use('plot_style.txt')
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.spines['bottom'].set_position(('data',0))
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_major_locator(MultipleLocator(tick_l))
    ax.xaxis.set_minor_locator(MultipleLocator(tick_s))
    ax.xaxis.set_major_formatter('{x:.0f}')
    if percent:
        ax.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1, decimals=None, symbol='%', is_latex=False))
    ax.tick_params(which='minor', width=1.0, length=3, labelsize=8)
    ax.axis([x1,y1,x2,y2])
    ax.set_xlabel(x_lab)
    ax.set_ylabel(y_lab)


############################################################
# Adjust column widths of Excel output file
def column_adjust(df, writer, sheet_name, max_width):
    worksheet = writer.sheets[sheet_name]
    #Iterate through each column and set the width == the max length in that column. A padding length of 2 is also added.
    for col_idx, col_name in enumerate(df.columns):
        col_name_width = len(col_name)
        # If column has names, identify length of longest name
        if isinstance(df.iloc[0][col_name], str):
            column_width = df[col_name].astype(str).str.len().max()
        else :
            column_width = 0
        # Determine if an item length is longer than column name length
        column_width = max(column_width, col_name_width) 
        if column_width > max_width:
            column_width = max_width
        if column_width < 6:
            column_width = 6
        worksheet.set_column(col_idx, col_idx, column_width)


############################################################
# Merge duplicate mutations by averaging ddG values

def remove_dups(df):
    pdb = df.pdb_id
    pdb_set = list(set(pdb))
    for pdbi in pdb_set:
        df_pdb = df[df.pdb_id==pdbi]
        pos = df_pdb.position
        pos_set = list(set(pos))
        for posi in pos_set:
            df_pdbp = df_pdb[df_pdb.position==posi]
            wild = df_pdbp.wild_type
            wild_set = list(set(wild))
            for wildi in wild_set:
                df_pdbp_w = df_pdbp[df_pdbp.wild_type==wildi]
                mut = df_pdbp_w.mutation
                mut_set = list(set(mut))
                for muti in mut_set:
                    df_pdbp_wm = df_pdbp_w[df_pdbp_w.mutation==muti]
                    ddg_av = df_pdbp_wm["ddG"].mean()
                    df.loc[(df["pdb_id"]==pdbi) & (df["position"]==posi) 
                           & (df["wild_type"]==wildi) & (df["mutation"]==muti), "ddG"] = ddg_av 
    df.drop_duplicates(keep='first', inplace=True)


#########################################################
# Add relavive asa values

def add_rel_asa(df):
    max_asa = pd.read_csv('Data Files/asa_max.csv')
    df["rel_asa"] = ""

    for ind in df.index:
        max_value = max_asa.loc[max_asa["wild_type"]==df["wild_type"][ind]]["max_asa"]
        df['rel_asa'][ind] = df['asa'][ind] / max_value.values[0]


#########################################################
# Add protein length

def add_len(df):
    df['len'] = ""

    for ind in df.index:
        prot_seq = df["sequence"][ind]
        prot_length = len(prot_seq)
        df['len'][ind] = prot_length


#########################################################
# Add destabilizing %, average ddG, and searchability
# Save as Excel sheet and print ddG histograms for both cutoffs

def destab_analysis(mut_data, sample, cutoff_1, cutoff_2, writer):
    # Determine destabilizing percentages for all data for both cutoffs 
    ddg = mut_data.ddG 
    stab = list(filter(lambda x: x < -cutoff_1, ddg))
    neut = list(filter(lambda x: x >= -cutoff_1, ddg))
    neut = list(filter(lambda x: x <= cutoff_1, neut))
    dest = list(filter(lambda x: x > cutoff_1, ddg))
    stab_per1 = len(stab)/len(ddg)
    neut_per1 = len(neut)/len(ddg)
    dest_per_tot1 = len(dest)/len(ddg)
    
    stab = list(filter(lambda x: x < -cutoff_2, ddg))
    neut = list(filter(lambda x: x >= -cutoff_2, ddg))
    neut = list(filter(lambda x: x <= cutoff_2, neut))
    dest = list(filter(lambda x: x > cutoff_2, ddg))
    stab_per2 = len(stab)/len(ddg)
    neut_per2 = len(neut)/len(ddg)
    dest_per_tot2 = len(dest)/len(ddg)

    # Calculate destabilizing percentage and average destabilization 
    #     for each protein and add to new dataframe
    # Determine searchability for each protein and add to new dataframe
    pid = mut_data.pdb_id
    pid_unique = list(set(pid))
    pid_sample = []
    pname = mut_data.protein_name
    pname_all = []
    pname_sample = []
    summary_data = pd.DataFrame()
    sample_data = pd.DataFrame()
    summary_data["pdb_id"] = pid_unique
    
    dest_all1 = []
    dest_all2 = []
    dest_sample1 = []
    dest_sample2 = []
    ddGav_all = []
    ddGav_sample = []
    len_all = []
    len_sample = []
    muts = []
    muts_sample = []
    searchable_prot1 = []
    searchable_prot2 = []
    searchable_meta1 = []
    searchable_meta2 = []
    for pidi in pid_unique:
        mutdata_id = mut_data[mut_data.pdb_id==pidi]
        pname = mutdata_id.iloc[0]['protein_name']
        pname_all.append(pname)
        plen = mutdata_id.iloc[0]['len']
        len_all.append(plen)
        ddg_id = mutdata_id.ddG
        
        dest = list(filter(lambda x: x > cutoff_1, ddg_id))
        dest_per1 = len(dest)/len(ddg_id)
        dest_all1.append(dest_per1)
        
        dest = list(filter(lambda x: x > cutoff_2, ddg_id))
        dest_per2 = len(dest)/len(ddg_id)
        dest_all2.append(dest_per2)
        
        ddGav =  mutdata_id["ddG"].mean()
        ddGav_all.append(ddGav)
        
        mut_num = len(ddg_id)
        muts.append(mut_num)
        if len(ddg_id) >= sample:
            pid_sample.append(pidi)
            pname_sample.append(pname)
            
            dest_sample1.append(dest_per1)
            dest_sample2.append(dest_per2)
            
            ddGav_sample.append(ddGav)
            len_sample.append(plen)
            muts_sample.append(mut_num)
            if dest_per1 == 0 :
                searchable_prot1.append(False)
                searchable_meta1.append(False)
            elif plen < -20/math.log((1-dest_per1), 10):
                searchable_prot1.append(True)
                searchable_meta1.append(True)
            elif plen < -40/math.log((1-dest_per1), 10):
                searchable_prot1.append(True)
                searchable_meta1.append(False)
            else :
                searchable_prot1.append(False)
                searchable_meta1.append(False)
                
            if dest_per2 == 0 :
                searchable_prot2.append(False)
                searchable_meta2.append(False)
            elif plen < -20/math.log((1-dest_per2), 10):
                searchable_prot2.append(True)
                searchable_meta2.append(True)
            elif plen < -40/math.log((1-dest_per2), 10):
                searchable_prot2.append(True)
                searchable_meta2.append(False)
            else :
                searchable_prot2.append(False)
                searchable_meta2.append(False)
                            
    cutoff_str1 = "({:.1f})".format(cutoff_1)
    des_per_id1 = "des_per  " + cutoff_str1 
    cutoff_str2 = "({:.1f})".format(cutoff_2)
    des_per_id2 = "des_per  " + cutoff_str2
    summary_data["protein_name"] = pname_all
    summary_data[des_per_id1] = dest_all1
    summary_data[des_per_id2] = dest_all2
    summary_data["ddG_av"] = ddGav_all
    summary_data["len"] = len_all
    summary_data["sample"] = muts
    
    searchable_id_p1 = "search prot  " + cutoff_str1 
    searchable_id_m1 = "search meta  " + cutoff_str1 
    searchable_id_p2 = "search prot  " + cutoff_str2
    searchable_id_m2 = "search meta  " + cutoff_str2
    sample_data["pub_id"] = pid_sample
    sample_data["protein_name"] = pname_sample
    sample_data[des_per_id1] = dest_sample1
    sample_data[des_per_id2] = dest_sample2
    sample_data["ddG_av"] = ddGav_sample
    sample_data["len"] = len_sample
    sample_data["sample"] = muts_sample
    sample_data[searchable_id_p1] = searchable_prot1
    sample_data[searchable_id_m1] = searchable_meta1
    sample_data[searchable_id_p2] = searchable_prot2
    sample_data[searchable_id_m2] = searchable_meta2
    searchable_per_p1 = searchable_prot1.count(True)/len(searchable_prot1)
    searchable_per_m1 = searchable_meta1.count(True)/len(searchable_meta1)
    searchable_per_p2 = searchable_prot2.count(True)/len(searchable_prot2)
    searchable_per_m2 = searchable_meta2.count(True)/len(searchable_meta2)

    # Write protein id, name, destabilizing %, average ddG, and searchable to to "mut_out [SAMPLE].xlsx
    summary_data.to_excel(writer, sheet_name='Proteins', index=False)
    workbook = writer.book
    worksheet = writer.sheets['Proteins']
    format1 = workbook.add_format({'num_format': '0%'})
    format2 = workbook.add_format({'num_format': '#,##0.00'})
    worksheet.set_column('C:D', None, format1)
    worksheet.set_column('E:E', None, format2)
    column_adjust(summary_data, writer, 'Proteins', 30)

    # Write sample data to "mut_out [SAMPLE].xlsx"
    sample_data.to_excel(writer, sheet_name='Samples', index=False)
    worksheet = writer.sheets['Samples']
    worksheet.set_column('C:D', None, format1)
    worksheet.set_column('E:E', None, format2)
    column_adjust(sample_data, writer, 'Samples', 30)
    
    # Add sheet with maximum searchable protein lengths for different percentiles of destabilizing %
    search_per = pd.DataFrame()
    percs = []
    perc_vals1 = []
    perc_vals2 = []
    max_len_p1 = []
    max_len_m1 = []
    max_len_p2 = []
    max_len_m2 = []
    dest_sample1.sort()
    dest_sample2.sort()
    for percentile in range(10, 100, 10):
        percs.append(percentile)
        dest_per1 = dest_sample1[math.floor(len(dest_sample1)*percentile*.01)]
        dest_per2 = dest_sample2[math.floor(len(dest_sample1)*percentile*.01)]
        perc_vals1.append(dest_per1)
        perc_vals2.append(dest_per2)
        max_len_p1.append(-40/math.log((1-dest_per1), 10))
        max_len_m1.append(-20/math.log((1-dest_per1), 10))
        max_len_p2.append(-40/math.log((1-dest_per2), 10))
        max_len_m2.append(-20/math.log((1-dest_per2), 10))
    dest_label1 = "dest %  " + cutoff_str1 
    dest_label2 = "dest %  " + cutoff_str2 
    perc_id_p1 = "max len prot  " + cutoff_str1 
    perc_id_m1 = "max len meta  " + cutoff_str1 
    perc_id_p2 = "max len prot  " + cutoff_str2
    perc_id_m2 = "max len meta  " + cutoff_str2
    search_per['percentile'] = percs
    search_per[dest_label1] = perc_vals1
    search_per[perc_id_p1] = max_len_p1
    search_per[perc_id_m1] = max_len_m1
    search_per[dest_label2] = perc_vals2
    search_per[perc_id_p2] = max_len_p2
    search_per[perc_id_m2] = max_len_m2
    
    search_per.to_excel(writer, sheet_name='Percentiles', index=False)  
    worksheet = writer.sheets['Percentiles']
    format3 = workbook.add_format({'num_format': '#,##0'})
    worksheet.set_column('B:B', None, format1)
    worksheet.set_column('C:D', None, format3)
    worksheet.set_column('E:E', None, format1)
    worksheet.set_column('F:G', None, format3)
    column_adjust(search_per, writer, 'Percentiles', 30)

    # Save coumulative data to "mut_out [SAMPLE].xlsx"
    stab_per_label1 = "Stabilizing %  " + cutoff_str1 
    neut_per_label1 = "Neutral %  " + cutoff_str1 
    dest_per_label1 = "Destabilizing %  " + cutoff_str1 
    searchable_label_p1 = "Protozoan Searchable %  " + cutoff_str1 
    searchable_label_m1 = "Metazoan Searchable %  " + cutoff_str1 
    stab_per_label2 = "Stabilizing %  " + cutoff_str2 
    neut_per_label2 = "Neutral %  " + cutoff_str2 
    dest_per_label2 = "Destabilizing %  " + cutoff_str2 
    searchable_label_p2 = "Protozoan Searchable %  " + cutoff_str2
    searchable_label_m2 = "Metazoan Searchable %  " + cutoff_str2 
    output_data = [[stab_per_label1, stab_per1], [neut_per_label1, neut_per1], 
                   [dest_per_label1, dest_per_tot1],
                   [searchable_label_p1, searchable_per_p1], 
                   [searchable_label_m1, searchable_per_m1], ["  ", "  "],
                   [stab_per_label2, stab_per2], [neut_per_label2, neut_per2], 
                   [dest_per_label2, dest_per_tot2],
                   [searchable_label_p2, searchable_per_p2], 
                   [searchable_label_m2, searchable_per_m2]]            
    df_data = pd.DataFrame(output_data, columns = ['Variable', 'Value'])
    
    df_data.to_excel(writer, sheet_name='Combined', index=False)
    worksheet = writer.sheets['Combined']
    worksheet.set_column('B:B', None, format1)
    column_adjust(df_data, writer, 'Combined', 25)
    
    # Add number of proteins in analysis to commulative data sheet
    worksheet.set_column(3, 3, 12)
    worksheet.write(1, 3, 'ddG Proteins')
    worksheet.write(1, 4, len(pid_sample))
    
    # Plot figure for histograms of destabilizing % for both cutoffs
    plot1 = plt.figure(1)
    fig, axs = plt.subplots(nrows=1, ncols=2)
    fig.set_size_inches(5.5, 3)
    plt.subplots_adjust(left=0.1,
                        bottom=0.15, 
                        right=0.95, 
                        top=0.97, 
                        wspace=0.2,     
                        hspace=0.1)
    
    # Distribution of destabilizing % for cutoff_1
    xaxis_title = "Destabilizing % ({:.1f} kcal/mol)".format(cutoff_1)
    apply_styles(axs[0], 0, 1.05, 0, 10, .2, .1, xaxis_title, 'Frequency', False)
    axs[0].xaxis.set_major_formatter(mtick.PercentFormatter(xmax=1, decimals=None, 
                                                           symbol='%', is_latex=False))
#    axs[0].text(0.9, 1.0, '(a)', horizontalalignment='left', fontsize=8,
#                verticalalignment='top', transform=axs[0].transAxes)
    axs[0].hist(dest_sample1, bins=10, range=(0,1), rwidth=0.8, density=False, align='mid')
    
    # Distribution of destabilizing % for cutoff_2
    xaxis_title = "Destabilizing % ({:.1f} kcal/mol)".format(cutoff_2)
    apply_styles(axs[1], 0, 1.05, 0, 10, .2, .1, xaxis_title, 'Frequency', False)
    axs[1].xaxis.set_major_formatter(mtick.PercentFormatter(xmax=1, decimals=None, 
                                                           symbol='%', is_latex=False))
#    axs[1].text(0.9, 1.0, '(b)', horizontalalignment='left', fontsize=8,
#                verticalalignment='top', transform=axs[1].transAxes)
    axs[1].hist(dest_sample2, bins=10, range=(0,1), rwidth=0.8, density=False, align='mid')
    
    # Add figure title
    fig.subplots_adjust(top=0.82)

    # Save figure of destabilizing % histograms
    file_name = "Figures/Histograms " + str(sample) + ".png"
    plt.savefig(file_name)
    
    # Plot figures of length versus destabilizing %
    plot2 = plt.figure(2)
    fig, axs = plt.subplots(nrows=1, ncols=2)
    fig.set_size_inches(5.5, 3)
    plt.subplots_adjust(left=0.15,
                        bottom=0.15, 
                        right=0.9, 
                        top=0.9, 
                        wspace=0.4,     
                        hspace=0.1)
    
    # Length versus destabilizing % for cutoff = cutoff_1
    apply_styles(axs[0], 0, 2700, 0, 1.05, 500, 250, 'Length', 'Destabilizing %', True)
    axs[0].text(0.05, 1.0, '(a)', horizontalalignment='left', fontsize=8,
                verticalalignment='top', transform=axs[0].transAxes)
    axs[0].scatter(len_sample, dest_sample1)
    
    # Length (expanded) versus destabilizing % for cutoff = cutoff_1 
    apply_styles(axs[1], 0, 500, 0, 1.05, 100, 50, 'Length (Contracted)', 'Destabilizing %', True)
    axs[1].text(0.05, 1.0, '(b)', horizontalalignment='left', fontsize=8,
                verticalalignment='top', transform=axs[1].transAxes)
    axs[1].scatter(len_sample, dest_sample2)
    
    # Add figure title
    fig.subplots_adjust(top=0.8)
    
    # Save figure of Relationships
    file_name = "Figures/Len-Dest " + str(sample) + " " + str(cutoff_1) + ".png"
    plt.savefig(file_name)


############################################################
# Calculate core % and average relative asa for each protein
# Plot relationships and save data to file

def relationship_analysis(mut_data, sample, cutoff_1, cutoff_2, writer):
    # Calculate core % and average relative asa for each protein
    sample_data = pd.DataFrame()
    pid = mut_data.pdb_id
    pid_unique = list(set(pid))
    pid_sample = []
    pname = mut_data.protein_name
    pname_sample = []
    
    dest_sample1 = []
    dest_sample2 = []
    ddGav_sample = []
    
    muts_sample = []
    core_sample = []
    len_sample = []
    rasa_ave = []
    for pidi in pid_unique:
        mutdata_id = mut_data[mut_data.pdb_id==pidi]
        ddg_id = mutdata_id.ddG
        rasa_id = mutdata_id.rel_asa
        if len(ddg_id) >= sample: 
            pid_sample.append(pidi)
            pname = mutdata_id.iloc[0]['protein_name']
            pname_sample.append(pname)
            
            dest = list(filter(lambda x: x > cutoff_1, ddg_id))
            dest_per = len(dest)/len(ddg_id)
            dest_sample1.append(dest_per)
            
            dest = list(filter(lambda x: x > cutoff_2, ddg_id))
            dest_per = len(dest)/len(ddg_id)
            dest_sample2.append(dest_per)
            
            ddGav =  mutdata_id["ddG"].mean()
            ddGav_sample.append(ddGav)
            
            core = list(filter(lambda x: x < .25, rasa_id))
            core_per = len(core)/len(rasa_id)
            core_sample.append(core_per)
            rasa_ave.append(st.mean(rasa_id))
            len_value = mutdata_id.iloc[0]['len']
            len_sample.append(len_value)
            muts_sample.append(len(ddg_id))
            
    cutoff_str1 = "({:.1f})".format(cutoff_1)
    des_per_id1 = "des_per  " + cutoff_str1 
    cutoff_str2 = "({:.1f})".format(cutoff_2)
    des_per_id2 = "des_per  " + cutoff_str2
    sample_data["pid"] = pid_sample
    sample_data["protein_name"] = pname_sample
    sample_data[des_per_id1] = dest_sample1
    sample_data[des_per_id2] = dest_sample2
    sample_data["ddG_av"] = ddGav_sample
    sample_data["rasa_av"] = rasa_ave
    sample_data["core %"] = core_sample
    sample_data["len"] = len_sample
    sample_data["sample"] = muts_sample

    
    # Plot relationships between destabilizing %, ASA, and length
    fig, axs = plt.subplots(nrows=1, ncols=2)
    fig.set_size_inches(5.5, 3)
    plt.subplots_adjust(left=0.15,
                        bottom=0.15, 
                        right=0.9, 
                        top=0.9, 
                        wspace=0.4,     
                        hspace=0.1)
    
    # Core % versus destabilizing % for cutoff = cutoff_1
    apply_styles(axs[0], 0, 1.05, 0, 1, .2, .1, 'Core %', 'Destabilizing %', True)
    axs[0].xaxis.set_major_formatter(mtick.PercentFormatter(xmax=1, decimals=None, 
                                                           symbol='%', is_latex=False))
    axs[0].text(0.05, 1.0, '(a)', horizontalalignment='left', fontsize=8,
                verticalalignment='top', transform=axs[0].transAxes)
    axs[0].scatter(core_sample, dest_sample1)
    
    # Length versus core %
    apply_styles(axs[1], 0, 500, 0, 1.05, 100, 50, 'Length', 'Core %', True)
    axs[1].text(0.05, 1.0, '(b)', horizontalalignment='left', fontsize=8,
                verticalalignment='top', transform=axs[1].transAxes)
    axs[1].scatter(len_sample, core_sample)
    
    # Add figure title
    fig.subplots_adjust(top=0.8)

      
    # Save figure of Relationships
    file_name = "Figures/Core-Len " + str(sample) + " " + str(cutoff_1) + ".png"
    plt.savefig(file_name)

    # Write protein id, name, destabilizing %, average ddG, average relative ASA, and core % to writer
    sample_data.to_excel(writer, sheet_name='ASA Samples', index=False)
    workbook = writer.book
    worksheet = writer.sheets['ASA Samples']
    format1 = workbook.add_format({'num_format': '0%'})
    format2 = workbook.add_format({'num_format': '#,##0.00'})
    worksheet.set_column('C:D', None, format1)
    worksheet.set_column('E:E', None, format2)
    worksheet.set_column('F:G', None, format1)
    column_adjust(sample_data, writer, 'ASA Samples', 30)
    
    # Add number of proteins in analysis to commulative data sheet
    worksheet = writer.sheets['Combined']
    worksheet.write(2, 3, 'ASA Proteins')
    worksheet.write(2, 4, len(pid_sample))


############################################################
# Create destabilizing distributions and destabilizing percentages for Tawfik data

def mut_dists(mut_data, cutoff_1, cutoff_2):
    # Create ddG distributions
    x_calc = []
    prot_dists = {}
    stab_per1 = []
    neut_per1 = []
    dest_per1 = []
    stab_per2 = []
    neut_per2 = []
    dest_per2 = []
    for prot_num in range(16):
        ddg_dist = []
        prot = mut_data.iloc[prot_num]
        dist_func = lambda x:gauss(x, prot['mean1'], prot['sig1'], prot['P1']) + gauss(x, prot['mean2'], prot['sig2'], 1-prot['P1'])
        for pos in range(81):
            x = -6 + pos*.4
            if prot_num==0:
                x_calc.append(x)
            dist_value = dist_func(x)
            ddg_dist.append(dist_value)
        prot_dists[prot_num] = ddg_dist
        stab_per1.append(scipy.integrate.quad(dist_func, -4, -cutoff_1)[0])
        neut_per1.append(scipy.integrate.quad(dist_func, -cutoff_1, cutoff_1)[0])
        dest_per1.append(scipy.integrate.quad(dist_func, cutoff_1, 10)[0])
        stab_per2.append(scipy.integrate.quad(dist_func, -4, -cutoff_2)[0])
        neut_per2.append(scipy.integrate.quad(dist_func, -cutoff_2, cutoff_2)[0])
        dest_per2.append(scipy.integrate.quad(dist_func, cutoff_2, 10)[0])

    
    # Add stabilizing, neutral, and destabilizing % to dataframe for both cutoffs
    mut_data["stab %  (0.5)"] = stab_per1
    mut_data["neut %  (0.5)"] = neut_per1
    mut_data["dest %  (0.5)"] = dest_per1
    mut_data["stab %  (1.0)"] = stab_per2
    mut_data["neut %  (1.0)"] = neut_per2
    mut_data["dest %  (1.0)"] = dest_per2
    
    # Plot ddG distributions
    plot1 = plt.figure(1)
    abb = mut_data.Abb
    fig, axs = plt.subplots(nrows=2, ncols=1)
    fig.set_size_inches(5, 7)
    plt.subplots_adjust(left=0.2,
                        bottom=0.07, 
                        right=0.9, 
                        top=0.95, 
                        wspace=0.4, 
                        hspace=0.3)
    
    markers = ['*', '+', 'x', '2', '4', 'P', 'p', 'X', '*', '+', 'x', '2', '4', 'P', 'p', 'X']
    apply_styles(axs[0], -6, 10, 0, .5, 2, 1, '$\Delta\Delta$G (kcal/mol)', 'Fraction (%)', True)
    for plotn in range(8):
        axs[0].plot(x_calc,prot_dists[plotn], marker=markers[plotn], label=abb[plotn])
    axs[0].legend()
   
    
    apply_styles(axs[1], -6, 10, 0, .5, 2, 1, '$\Delta\Delta$G (kcal/mol)', 'Fraction (%)', True)
    for plotn in range(8,16):
        axs[1].plot(x_calc,prot_dists[plotn], marker=markers[plotn],label=abb[plotn])
    axs[1].legend()
  

    
    # Save figure of Tawfik protein ddG distribution functions
    file_name = "Figures/Tawfik Proteins.png"
    plt.savefig(file_name)
    