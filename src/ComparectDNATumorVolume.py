'''
Author: Juber Patel
'''

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import date
from matplotlib import rcParams
from matplotlib.ticker import ScalarFormatter
from pandas.api.types import CategoricalDtype
import sys
import colorcet as cc
import scipy.stats as stats
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter




def main():
    
    d = "/Users/patelj1/current/Beryl-17-180/analysis-new/tumor-volume-correlation/"

    os.chdir(d)

    comparisonFile = "tumor-volume-ctDNA-fraction-comparison.tsv"
    
    df = pd.read_csv(comparisonFile, sep="\t", index_col=False, keep_default_na=False, low_memory=False)

    '''
    # compute the slopes
    for p in df["Patient"].unique():
        
        df.loc[df["Patient"] == p, "ctDNAFractionDiff"] = df[(df["Patient"] == p) & (df["Timepoint"]== "TP2")]["ctDNAFraction"].values[0] - df[(df["Patient"] == p) & (df["Timepoint"]== "Base")]["ctDNAFraction"].values[0]
        
        df.loc[df["Patient"] == p, "ctDNAFractionInterval"] = df[(df["Patient"] == p) & (df["Timepoint"]== "TP2")]["cfDNA Days"].values[0] - df[(df["Patient"] == p) & (df["Timepoint"]== "Base")]["cfDNA Days"].values[0]

        df.loc[df["Patient"] == p, "VolumeDiff"] = df[(df["Patient"] == p) & (df["Timepoint"]== "TP2")]["Volume (cm3)"].values[0] - df[(df["Patient"] == p) & (df["Timepoint"]== "Base")]["Volume (cm3)"].values[0]

        df.loc[df["Patient"] == p, "VolumeInterval"] = df[(df["Patient"] == p) & (df["Timepoint"]== "TP2")]["CT Scan Days"].values[0] - df[(df["Patient"] == p) & (df["Timepoint"]== "Base")]["CT Scan Days"].values[0]

    
    df["ctDNASlope"] = df["ctDNAFractionDiff"] / df["ctDNAFractionInterval"]
    df["VolumeSlope"] = df["VolumeDiff"] / df["VolumeInterval"]
    '''

    # compute fold change

    # change this!!
    df = df[(df["ctDNAFraction"]!=0) & (df["Volume (cm3)"]!=0)]
    
    for p in df["Patient"].unique():
        
        df.loc[df["Patient"] == p, "ctDNAFractionFC"] = df[(df["Patient"] == p) & (df["Timepoint"]== "TP2")]["ctDNAFraction"].values[0]/df[(df["Patient"] == p) & (df["Timepoint"]== "Base")]["ctDNAFraction"].values[0]

        df.loc[df["Patient"] == p, "VolumeFC"] = df[(df["Patient"] == p) & (df["Timepoint"]== "TP2")]["Volume (cm3)"].values[0]/df[(df["Patient"] == p) & (df["Timepoint"]== "Base")]["Volume (cm3)"].values[0]



    df.to_csv("slope-fc.tsv", sep="\t", index=False)

    df.drop_duplicates(subset=["Patient"], inplace=True)

    # Pearson correlation
    print(stats.pearsonr(df["ctDNAFractionFC"], df["VolumeFC"]))

    #stats.wilcoxon

    # compute correlation coefficient


    # plot the slopes

    sns.set_theme(style="whitegrid")

    fig, ax = plt.subplots(figsize=(10, 10))
    
    plt.figure(dpi=300)
    plt.tight_layout()
    
    g = sns.scatterplot(data = df, x="VolumeFC", y="ctDNAFractionFC")
    g.set_xlabel(xlabel="Base:8wk tumor volume fold change")
    g.set_ylabel(ylabel="Base:2wk ctDNA fraction fold change")

    plt.tight_layout()
    plt.savefig("tumor-volume-ctDNA-fraction-scatterplot.pdf")
    plt.close()


    # plot the fc
    # convert to long form
    df = df.melt(id_vars=['Patient'], value_vars=['ctDNAFractionFC', 'VolumeFC'], var_name='Method', value_name='Fold Change')

    df.to_csv("t.tsv", sep="\t", index=False)
    

    fig, ax = plt.subplots(figsize=(10, 10))

    plt.figure(dpi=300)
    plt.tight_layout()
    
    g = sns.barplot(data = df, x="Patient", y="Fold Change", hue="Method")
    g.axhline(1.0, color='black')
    g.set_yscale('log', base=2)
    #ax.yaxis.set_major_formatter(FormatStrFormatter('%.4f'))
    plt.xticks(rotation=90)
    

    plt.tight_layout()
    plt.savefig("fc.pdf")
    plt.close()









if __name__ == "__main__":
    main()
    
