'''
Created on March 30, 2023

@author: Juber Patel
'''


import os
import sys
import time
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats as stats




'''
Process data from all cells for all mutations.
This is where filtering happens and cell fractions for alleles are computed.
'''
def main():
    

    d = "/Users/patelj1/current/Autopsy/analysis/filtered-new/MSK-AB-0002/MissionBio/analysis"

    os.chdir(d)

    genotypingFile = "genotypes-selected-columns.tsv"
    bulkCCFsFile = "MSK-AB-0002-absolute-filtered-annotated-ccfs.txt"

    truncalDriverMutations = ["CDH1.p.W156Tfs*60;LikelyOncogenic;;16:68842402:G:GAC;;;Cancer", "TP53.p.N239D;LikelyOncogenic;HOTSPOT;17:7577566:T:C;;;Cancer"]

    sampleName = "T23_1-left-top-back-skin-met-1"


    genotypes = pd.read_csv(genotypingFile, sep="\t", index_col=False, keep_default_na=False, low_memory=False)

    min_depth = 15
    min_AF_presence = 0.10

    genotypes = genotypes[['Tumor_Sample_Barcode', 'MYID', 'Waltz_total_t_depth', 'Waltz_total_t_alt_count']]
                          

    genotypes = genotypes[genotypes['Waltz_total_t_depth'] >= min_depth]
    genotypes["AF"] = genotypes["Waltz_total_t_alt_count"] / genotypes["Waltz_total_t_depth"]

    tumorCells = genotypes[(genotypes["MYID"].isin(truncalDriverMutations)) & (genotypes["AF"] >= min_AF_presence)]["Tumor_Sample_Barcode"].unique().tolist()

    tumorCellsCount = len(tumorCells)
    print(f"Found {tumorCellsCount} cells with truncal driver mutations")

    # only keep tumor cells
    genotypes = genotypes[genotypes["Tumor_Sample_Barcode"].isin(tumorCells)]

    genotypes["Present"] = ((genotypes["AF"] >= min_AF_presence) | ((genotypes["AF"]>=0.05) & (genotypes["Waltz_total_t_alt_count"]>=10)))

    genotypes["PresentIn"] = genotypes.groupby("MYID")["Present"].transform("sum")
    genotypes["CCF_sc"] = genotypes["PresentIn"] / tumorCellsCount

    genotypes.to_csv("sc-muts.tsv", sep="\t", index=False)

    ccfs_sc = genotypes.drop_duplicates(subset=["MYID"])
    ccfs_sc = ccfs_sc[["MYID", "CCF_sc"]]
    


    # get bulk CCFs
    bulkCCFs = pd.read_csv(bulkCCFsFile, sep="\t", index_col=False, keep_default_na=False, low_memory=False)

    bulkCCFs = bulkCCFs[bulkCCFs["Tumor_Sample_Barcode"] == sampleName]

    ccfs = pd.melt(bulkCCFs, id_vars="Tumor_Sample_Barcode", var_name='MYID', value_name='CCF_bulk')
    ccfs = ccfs.merge(ccfs_sc, on="MYID", how="left")
    ccfs = ccfs[~ccfs["CCF_sc"].isna()]
    ccfs["CCF_sc"] = ccfs["CCF_sc"].astype(float)
    ccfs["CCF_bulk"] = ccfs["CCF_bulk"].astype(float)

    ccfs.to_csv("ccfs.tsv", sep="\t", index=False)



    # plot the CCFs

    # plot histogram of position error
    sns.set_style("whitegrid")
    
    
    # plot per-cell AF for truncal driver mutations
    fig = plt.figure(dpi=300)
    plt.tight_layout()
    t = genotypes[genotypes["MYID"].isin(truncalDriverMutations)].copy()
    g = sns.displot(data=t, x="AF", hue="MYID", kind="kde")
    sns.move_legend(g, "center")
    plt.title("Distribution of in-cell AFs for truncal driver mutations")
    plt.tight_layout()
    plt.savefig("truncal-driver-per-cell-AF-kde.pdf")

    fig = plt.figure(dpi=300)
    plt.tight_layout()
    g = sns.displot(data=t, x="AF", bins=20, hue="MYID", multiple="dodge")
    sns.move_legend(g, "center")
    plt.title("Distribution of in-cell AFs for truncal driver mutations")
    plt.tight_layout()
    plt.savefig("truncal-driver-per-cell-AF-hist.pdf")

    


    # plot bulk vs sc CCFs scatterplot
    print(ccfs.shape)
    # Pearson correlation
    print(stats.pearsonr(ccfs["CCF_bulk"], ccfs["CCF_sc"]))


    fig = plt.figure(dpi=300)
    plt.tight_layout()
    g = sns.scatterplot(data=ccfs, x="CCF_bulk", y="CCF_sc", alpha=0.5)
    #g.set_ylabel("Cumulative fraction of positions")
    #g.set_xlabel("Position error (fraction of cells with alt alleles)")
    plt.title(sampleName)
    plt.tight_layout()
    plt.savefig(sampleName + "-bulk-sc-ccfs.pdf")


        



    






    




if __name__ == '__main__':
    main()


