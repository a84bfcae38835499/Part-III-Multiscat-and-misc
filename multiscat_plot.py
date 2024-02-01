#%%
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 10 16:42:27 2023

This is a suggested method for plotting the output results of multiscat.

@author: SamLambrick
"""

import numpy as np
import pandas as pd
import seaborn as sns
import datetime
from matplotlib import pyplot as plt
# Default theme
sns.set_theme()

def import_multiscat(fname):
    """Impot standrad multiscat output into a pandas data frame."""
    d = pd.read_csv(fname, skiprows=7, delim_whitespace=True, 
                    header=None, names=['#','n1','n2','I'])
    d.drop(columns=['#'], inplace=True)
    return(d)

def calculate_entropy(dataf):
    H = 0
    vals = dataf.values
    for p in np.nditer(vals):
        if(np.isnan(p)):
            print("NaN found! skipping...")
        elif(p == 0):
            print("Zero found! skpping...")
        else:
            print("p = " + str(p))
            H += p * np.log(p)
    H /= -np.log(dataf.size)
    return(H)

d = import_multiscat('diffrac10001.out')


d2 = d.pivot(index='n1', columns='n2', values='I')


H = calculate_entropy(d2)
print("Entropy = " + str(H))

ax = sns.heatmap(d2, cmap='viridis', cbar_kws={'label' : '$P(n_1,n_2)$'})
ax.set_aspect('equal')
ax.set_xlabel('$n_1$')
ax.set_ylabel('$n_2$')

savestr = "Figures/Diffraction/" + datetime.datetime.now().strftime('Diffraction_%Y-%m-%d_%H-%M') + ".png"
plt.savefig(fname=savestr)
# %%

