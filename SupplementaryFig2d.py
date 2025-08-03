#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  1 09:38:14 2025

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

R0perc = Find("/path/to/tideR0_percentages.csv")
R0perc= R0perc.dropna() #remove empty rows
R0pval = Find("/path/to/tideR0_pvals.csv")
R0pval= R0pval.dropna() #remove empty rows

#where p>0.001 (non-significant indel) remove from consideration
R0perc.loc[:, ~R0perc.columns.isin(['replicate', 'indels'])] = R0perc.loc[:, ~R0perc.columns.isin(['replicate', 'indels'])].mask(R0pval > 0.001, 0)

#formating for plotting
melted_df = R0perc.melt(id_vars=['replicate', 'indels'],
                               var_name='condition',
                               value_vars= ['sgLAD_6hr',
                                      'sgLAD_24hr', 'sgLAD_48hr'],
                               value_name='percentage')
melted_df['percentage'] = melted_df['percentage'].astype(int)

#%%
#assign indel type 
def classify_indel(indel):
    if indel == 0:
        return "uncut"
    elif -4 <= indel <= 2:
        return "NHEJ"
    elif -20 <= indel <= -3:
        return "MMEJ"
    else:
        return None


melted_df['class'] = melted_df['indels'].apply(classify_indel)
melted_df = melted_df[melted_df['class'].notna()]

#find total % of sequences that are either NHEJ or MMEJ. (uncut/wt sequences are only ever "0")
summed_df = (
    melted_df
    .groupby(['condition', 'replicate', 'class'], as_index=False)
    .agg(summed_percentage=('percentage', 'sum'))
)


subset = summed_df.copy()

#Calculate total NHEJ + uncut per condition/replicate
total = subset.groupby(['condition', 'replicate'])['summed_percentage'].transform('sum')
#Calculate percentage within that total
subset['relative_percentage'] = (subset['summed_percentage'] / total) * 100
# Split 'condition' into two new columns: 'condition' and 'time'
subset[['condition', 'time']] = subset['condition'].str.split('_', n=1, expand=True)



#%%

hue_order = ['6hr', '24hr', '48hr']

custom_palette = {
    'MMEJ': '#D33873',
    'NHEJ': '#575757',
    'uncut': '#EFB54F'
}

sns.set_context("paper", font_scale=2.5)


g = sns.catplot(
    data=subset,
    x="class",
    y="relative_percentage",
    hue="class",
    col="time",
    kind="bar",
    palette=custom_palette,
    col_order=hue_order,
    errorbar= None,
    height=4,
    aspect=1.8,
    dodge=False,
    legend=False,
    alpha = 0.9
)

for ax, cond in zip(g.axes.flat, hue_order):
    data = subset[subset['time'] == cond]
    sns.stripplot(
        data=data,
        x="class",
        y="relative_percentage",
        ax=ax,
        palette=["white"],  
        alpha=0.8,
        size=14,
        edgecolor='black',
        linewidth=0.5,
        dodge=False,
        order=["MMEJ", "NHEJ", "uncut"]
    )

g.set_axis_labels("", "% of sequences")
g.set_titles("{col_name}")
plt.tight_layout()
plt.show()

plt.savefig("bypathwayrelative.svg", dpi=300, bbox_inches='tight')



