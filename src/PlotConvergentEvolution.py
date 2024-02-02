'''
Created on June 20, 2023

@author: Juber Patel
'''

import pandas as pd
import os
import seaborn as sb
import copy
import matplotlib.pyplot as plt
import sys
import seaborn as sns
from matplotlib import rcParams
from matplotlib.ticker import MaxNLocator



def main():
    
    d = "/Users/patelj1/current/Autopsy/analysis/filtered-new/convergent-evolution"

    os.chdir(d)

    muts = pd.read_csv("muts1.tsv", sep="\t", index_col=False, keep_default_na=False, low_memory=False)

    patients = muts["Patient"].unique()


    # make plots
    rcParams['figure.figsize'] = 5,40 
    sns.set_theme(style="whitegrid")
    plt.figure(dpi=300)
    plt.tight_layout()

    os.chdir("plots")

    for patient in patients:

        #if patient != "MSK-AB-0001":
        #    break

        print(patient)
        df = muts[(muts["Patient"]==patient)].copy()
        df.drop_duplicates(subset=["Gene"], inplace=True)
        df["MutsPer100Codons"] = df["MutsPer100Codons"].round(2)
        df = df[df["MutsPer100Codons"]>=0.5]
        
        print(df.shape[0])
        if df.shape[0] == 0:
            continue

        height = df.shape[0]/5
        if height < 2:
            height = 2

        rcParams['figure.figsize'] = 5,height

        df["GeneName"] = df["Gene"] + "-" + df["MutsPer100Codons"].astype(str)

        #df.to_csv("t.tsv", sep="\t", index=False)
        
        
        #p = sns.color_palette(cc.glasbey, n_colors=numPatients)
        title = patient + " Genes with highest mutations per 100 codons"
        g = sns.barplot(data = df, y="GeneName", x="GeneMutationsDistinct", color="blue")
        g.xaxis.set_label_position('top')
        g.xaxis.set_ticks_position("top")
        plt.title(title, x=0.1)

        #plt.xticks(rotation=90)
        ax = plt.gca()
        #ax.yaxis.set_major_formatter(FormatStrFormatter('%.4f'))

        #plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
        plt.tight_layout()
        plt.savefig(title + ".pdf")
        plt.close()
    









if __name__ == '__main__':
    main()


