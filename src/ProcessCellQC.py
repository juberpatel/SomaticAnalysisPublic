'''
Created on March 30, 2023

@author: Juber Patel
'''


import os
import sys
import time
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt





'''
process the cell QC file
'''
def main():
    
    d = "/Users/patelj1/current/SingleCellAnalysis/analysis"

    os.chdir(d)

    WBCCellQCFile = "WBC-control-A/cell-qc.tsv"
    ISHICellQCFile = "ISHI-HEC6/cell-qc.tsv"
    EEC131CellQCFile = "EEC131/cell-qc.tsv"

    cellQC = pd.read_csv(WBCCellQCFile, sep="\t", index_col=False, keep_default_na=False, low_memory=False)
    cellQC["Sample"] = "WBC-control-A (7114)"

    ishiCellQC = pd.read_csv(ISHICellQCFile, sep="\t", index_col=False, keep_default_na=False, low_memory=False)
    ishiCellQC["Sample"] = "ISHI-HEC6 (1820)"

    eec131CellQC = pd.read_csv(EEC131CellQCFile, sep="\t", index_col=False, keep_default_na=False, low_memory=False)
    eec131CellQC["Sample"] = "EEC131 (9702)"

    cellQC = pd.concat([cellQC, ishiCellQC, eec131CellQC], axis=0)

    print(len(cellQC))


    # plot noise
    fig = plt.figure(dpi=300)
    #fig.set_size_inches(25,60)
    #plt.rcParams.update({'font.size': 22})
    plt.tight_layout()

    g = sns.boxplot(data=cellQC, y="Noise", x="Sample")
    #plt.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)
    plt.tight_layout()
    plt.savefig("Noise.pdf")

    # plot coverage
    fig = plt.figure(dpi=300)
    plt.tight_layout()

    g = sns.boxplot(data=cellQC, y="Coverage", x="Sample")
    plt.tight_layout()
    plt.savefig("Coverage.pdf")

    # plot good positions
    fig = plt.figure(dpi=300)
    plt.tight_layout()

    g = sns.boxplot(data=cellQC, y="GoodPositions", x="Sample")
    plt.tight_layout()
    plt.savefig("GoodPositions.pdf")

    # plot coverage vs noise scatter plot
    fig = plt.figure(dpi=300)
    plt.tight_layout()

    g = sns.scatterplot(data=cellQC[cellQC["Sample"]=="WBC-control-A (7114)"], x="Coverage", y="Noise", color="grey", alpha=0.3, s=5)
    plt.tight_layout()
    plt.savefig("coverage-noise.pdf")

    # plot coverage vs good positions scatter plot
    fig = plt.figure(dpi=300)
    plt.tight_layout()

    g = sns.scatterplot(data=cellQC[cellQC["Sample"]=="WBC-control-A (7114)"], x="Coverage", y="GoodPositions", color="grey", alpha=0.3, s=5)
    plt.tight_layout()
    plt.savefig("coverage-good-positions.pdf")




    








if __name__ == '__main__':
    main()

