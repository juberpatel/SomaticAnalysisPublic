'''
Created on Nov 10, 2022

@author: Juber Patel
'''

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import date
from matplotlib import rcParams
from matplotlib.ticker import ScalarFormatter
from pandas.api.types import CategoricalDtype
import sys
import colorcet as cc




def main():
    
    d = "/Users/patelj1/current/RNAMediatedDNARepair/IMPACT"
    
    os.chdir(d)
    
    deletionsFile = "impact-all-solid-tumor-deletions-categorized.txt"
    
    
    df = pd.read_csv(deletionsFile, sep="\t", index_col=False, keep_default_na=False, low_memory=False)
    
    df = df[df["Mutation_Status"]=="SOMATIC"]
    
    df = df[df["t_alt_count"]!=""]
    df = df[df["t_ref_count"]!=""]
    
    df["t_ref_count"] = df["t_ref_count"].astype(int)
    df["t_alt_count"] = df["t_alt_count"].astype(int)
    
    df = df[df["t_ref_count"]>=0]
    
    df["DP"] = df["t_ref_count"] + df["t_alt_count"]
    df["AF"] = df["t_alt_count"]/df["DP"]
    
    df = df[df["DP"]>=100]
    df = df[df["DP"]<3000]
    
    df["Deletion_Length"] = (df["End_Position"] - df["Start_Position"]) + 1
         
    sns.set_theme(style="whitegrid")
    
    plt.figure(dpi=300)
    plt.tight_layout()
    
    #p = sns.color_palette(cc.glasbey, n_colors=len(df["Deletion_Category"].unique()))
    
    ax = sns.scatterplot(data=df[df["Deletion_Category"]!="intron deletion"], x="DP", y="AF", color="grey", alpha=0.1)
    g = sns.scatterplot(data=df[df["Deletion_Category"]=="intron deletion"], x="DP", y="AF", color="red",alpha=1.0,
                         ax=ax)
    
    plt.tight_layout()
    plt.savefig("impact-deletions-dp-af.pdf")
    plt.close()
    

    plt.figure(dpi=300)
    plt.tight_layout()
    
    
    ax = sns.scatterplot(data=df[df["Deletion_Category"]!="intron deletion"], x="Deletion_Length", y="AF", color="grey", 
                         alpha=0.1)
    g = sns.scatterplot(data=df[df["Deletion_Category"]=="intron deletion"], x="Deletion_Length", y="AF", color="red",
                        alpha=1.0, ax=ax)


    plt.tight_layout()
    plt.savefig("impact-deletions-length-af.pdf")
    plt.close()
    
    
    plt.figure(dpi=300)
    plt.tight_layout()

    sns.histplot(data=df[df["Deletion_Category"]!="intron deletion"], x="AF")

    plt.tight_layout()
    plt.savefig("others-af-histogram.pdf")
    plt.close()
    
    
    
    plt.figure(dpi=300)
    plt.tight_layout()

    sns.histplot(data=df[df["Deletion_Category"]=="intron deletion"], x="AF")

    plt.tight_layout()
    plt.savefig("wid-af-histogram.pdf")
    plt.close()
    
    
    df.to_csv("df.tsv", sep="\t", index=False)
    
    
    
    
    
    
    
    

if __name__ == '__main__':
    main()
