#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  1 10:02:52 2025

@author: joannafernandez
"""


import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os 

#import data file and name cols 
def Find(file):
    genes = pd.read_csv(file, delimiter=",")
    genes.rename(columns={'Unnamed: 0': 'replicate'}, inplace=True)
    genes.rename(columns={'Unnamed: 1': 'indels'}, inplace=True)

    return genes

siperc = Find("/path/to/siH2AX_percentages.csv")
siperc= siperc.dropna() #remove empty rows
sipval = Find("/path/to/siH2AX_pvals.csv")
sipval= sipval.dropna() #remove empty rows

#where p>0.001 (non-significant indel) remove from consideration
siperc.loc[:, ~siperc.columns.isin(['replicate', 'indels'])] = siperc.loc[:, ~siperc.columns.isin(['replicate', 'indels'])].mask(sipval > 0.001, 0)

#formating for plotting
singles = siperc.melt(id_vars=['replicate', 'indels'],
                               var_name='condition',
                               value_vars= ['siSCR_6hr', 'siSCR_12hr', 'siSCR_24hr',
                                      'siSCR_48hr', 'siH2AX_6hr', 'siH2AX_12hr', 'siH2AX_24hr',
                                      'siH2AX_48hr'],
                               value_name='percentage')
singles['percentage'] = singles['percentage'].astype(int)

#assign indel type 
def classify_indel(indel):
    if indel == 0:
        return "uncut"
    elif -4 <= indel <= 2:
        return "NHEJ"
    else:
        return None


singles['class'] = singles['indels'].apply(classify_indel)
singles = singles[singles['class'].notna()]

#find total % of sequences that are either NHEJ or MMEJ. (uncut/wt sequences are only ever "0")
summed_df = (
    singles
    .groupby(['condition', 'replicate', 'class'], as_index=False)
    .agg(summed_percentage=('percentage', 'sum'))
)


subset = summed_df[summed_df['class'].isin(['NHEJ', 'uncut'])].copy()

#Calculate total NHEJ + uncut per condition/replicate
total = subset.groupby(['condition', 'replicate'])['summed_percentage'].transform('sum')
#Calculate percentage within that total
subset['relative_percentage'] = (subset['summed_percentage'] / total) * 100
# Split 'condition' into two new columns: 'condition' and 'time'
subset[['condition', 'time']] = subset['condition'].str.split('_', n=1, expand=True)


sns.set_context("paper", font_scale=1)
hue_order = ['6hr', '12hr', '24hr', '48hr']
custom_palette = {
    'NHEJ': '#575757',
    'uncut': '#EFB54F'
}


time_order = ['6hr', '12hr', '24hr', '48hr']
pathway_order = ['siSCR', 'siH2AX']

yticks = [0, 20, 40, 60, 80, 100]

g = sns.catplot(
    data=subset,
    x="class",
    y="relative_percentage",
    hue="class",
    col="time",  
    row="condition",
    kind="bar",
    palette=custom_palette,
    col_order=hue_order,
    row_order=pathway_order,
    errorbar=None,
    height=4,
    aspect=1.8,
    dodge=False,
    legend=False,
    alpha=0.9
)

for (row_val, col_val), ax in g.axes_dict.items():
    subsety = subset[(subset['condition'] == row_val) & (subset['time'] == col_val)]

    # Stripplot overlay
    sns.stripplot(
        data=subsety,
        x="class",
        y="relative_percentage",
        hue="class",
        dodge=False,
        palette=["white"],
        size=10,
        edgecolor="black",
        linewidth=0.6,
        alpha=0.8,
        ax=ax
    )
    
    ax.axhline(y=50, color='black', linestyle='--', linewidth=1)
    
for ax in g.axes.flat:
    ax.yaxis.set_ticks_position('left') 
    ax.tick_params(axis='y', which='both', labelleft=True)  
    ax.set_yticks([0, 20, 40, 60, 80, 100]) 
    
    ax.xaxis.set_ticks_position('bottom')
    ax.tick_params(axis='x', which='both', labelbottom=True)

    # Optional: Draw a horizontal reference line at y=50
    ax.axhline(50, linestyle='--', color='grey', linewidth=1)

# Tight layout to avoid clipping
sns.despine()
plt.tight_layout()
plt.show()    
    

plt.savefig("010825h2ax.svg", dpi=300, bbox_inches='tight')  # Or .pdf/.png
#%%%

#now for ATMI
atmiperc = Find("/path/to/ATMi_percentages.csv")
atmiperc= atmiperc.dropna() #remove empty rows
atmipval = Find("/path/to/ATMi_pvals.csv")
atmipval= atmipval.dropna() #remove empty rows

#where p>0.001 (non-significant indel) remove from consideration
atmiperc.loc[:, ~atmiperc.columns.isin(['replicate', 'indels'])] = atmiperc.loc[:, ~atmiperc.columns.isin(['replicate', 'indels'])].mask(atmipval > 0.001, 0)

#formating for plotting
inhibitor = atmiperc.melt(id_vars=['replicate', 'indels'],
                               var_name='condition',
                               value_vars= ['DMSO_6hr', 'DMSO_12hr', 'DMSO_24hr',
                                      'DMSO_48hr', 'ATMi_6hr', 'ATMi_12hr', 'ATMi_24hr', 'ATMi_48hr'],
                               value_name='percentage')
inhibitor['percentage'] = inhibitor['percentage'].astype(int)

#assign indel type 
inhibitor['class'] = inhibitor['indels'].apply(classify_indel)
inhibitor = inhibitor[inhibitor['class'].notna()]

#find total % of sequences that are either NHEJ or MMEJ. (uncut/wt sequences are only ever "0")
summed_df = (
    inhibitor
    .groupby(['condition', 'replicate', 'class'], as_index=False)
    .agg(summed_percentage=('percentage', 'sum'))
)

subset = summed_df[summed_df['class'].isin(['NHEJ', 'uncut'])].copy()

#Calculate total NHEJ + uncut per condition/replicate
total = subset.groupby(['condition', 'replicate'])['summed_percentage'].transform('sum')
#Calculate percentage within that total
subset['relative_percentage'] = (subset['summed_percentage'] / total) * 100
# Split 'condition' into two new columns: 'condition' and 'time'
subset[['condition', 'time']] = subset['condition'].str.split('_', n=1, expand=True)


#format figure
sns.set_context("paper", font_scale=1)
hue_order = ['6hr', '12hr', '24hr', '48hr']
custom_palette = {
    'NHEJ': '#575757',
    'uncut': '#EFB54F'
}


time_order = ['6hr', '12hr', '24hr', '48hr']
pathway_order = ['DMSO', 'ATMi']  
yticks = [0, 20, 40, 60, 80, 100]

# Your existing catplot
g = sns.catplot(
    data=subset,
    x="class",
    y="relative_percentage",
    hue="class",
    col="time",  
    row="condition",
    kind="bar",
    palette=custom_palette,
    col_order=hue_order,
    row_order=pathway_order,
    errorbar=None,
    height=4,
    aspect=1.8,
    dodge=False,
    legend=False,
    alpha=0.9
)

# Overlay stripplot and horizontal line
for (row_val, col_val), ax in g.axes_dict.items():
    subsety = subset[(subset['condition'] == row_val) & (subset['time'] == col_val)]

    # Stripplot overlay
    sns.stripplot(
        data=subsety,
        x="class",
        y="relative_percentage",
        hue="class",
        dodge=False,
        palette=["white"],
        size=10,
        edgecolor="black",
        linewidth=0.6,
        alpha=0.8,
        ax=ax
    )
    
    ax.axhline(y=50, color='black', linestyle='--', linewidth=1)
    
for ax in g.axes.flat:
    ax.yaxis.set_ticks_position('left')
    ax.tick_params(axis='y', which='both', labelleft=True)  
    ax.set_yticks([0, 20, 40, 60, 80, 100])
    
    ax.xaxis.set_ticks_position('bottom')
    ax.tick_params(axis='x', which='both', labelbottom=True)

    # Optional: Draw a horizontal reference line at y=50
    ax.axhline(50, linestyle='--', color='grey', linewidth=1)

# Tight layout to avoid clipping
sns.despine()
plt.tight_layout()
plt.show()    
    

plt.savefig("010825atmi.svg", dpi=300, bbox_inches='tight')  # Or .pdf/.png



