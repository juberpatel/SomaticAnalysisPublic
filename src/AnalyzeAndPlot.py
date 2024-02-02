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

    genotypingFile = "genotypes-panel-snps.maf"
    bulkCCFsFile = "MSK-AB-0002-absolute-filtered-annotated-ccfs.txt"

    truncalDriverMutations = ["hotspot:CDH1.p.W156Tfs*60;LikelyOncogenic;;16:68842402:G:GAC;;;Cancer", "hotspot:TP53.p.N239D;LikelyOncogenic;HOTSPOT;17:7577566:T:C;;;Cancer"]

    sampleName = "T23_1-left-top-back-skin-met-1"

    min_depth = 15
    min_AF_presence = 0.025
    max_AF_absence = 0.02

    genotypes = pd.read_csv(genotypingFile, sep="\t", index_col=False, keep_default_na=False, low_memory=False)


    os.chdir("plots")

    sns.set_style("whitegrid")

    # overall coverage
    genotypes["avg_depth"] = genotypes.groupby("MYID")["Waltz_total_t_depth"].transform("mean")

    # remove the ridiculous TP53 mutations
    genotypes = genotypes[~genotypes["MYID"].str.contains("TCCCAAGACTTAGTACCTGAAGGGTGAAATATTCTCCATCCAGTGGTTTCTTCTTTGGCTGGGGAGAGGAGCTGGTGTTGTTGGGCAGTGCTAGGAAAGAGGCAAGGAAAGGTGATAAAAGTGAATCTGAGGCATAACTGCACCCTTGGTCTCCTCCACCGCTTCTTGTCCTGCTTGCTTACCTCGCTTAGTGCTCCCTGGGGGCAGCTCGTGGTGAGGC")]

    genotypes = genotypes[~genotypes["MYID"].str.contains("TAGGGCCAGGAAGGGGCTGAGGTCACTCACCTGGAGTGAGCCCTGCTCCCCCCTGGCTCCTT")]

    


    '''
    # plot coverage
    df = genotypes.drop_duplicates(subset=["MYID"]).copy()
    
    fig = plt.figure(dpi=300)
    plt.tight_layout()
    g = sns.ecdfplot(data=df, x="avg_depth")
    plt.title("Panel Coverage")
    plt.tight_layout()
    plt.savefig("panel-coverage.pdf")
    plt.close()

    fig = plt.figure(dpi=300)
    plt.tight_layout()
    g = sns.ecdfplot(data=df, x="avg_depth")
    g.set_xlim(0, 100)
    plt.title("Panel Coverage")
    plt.tight_layout()
    plt.savefig("panel-coverage-zoomed.pdf")
    plt.close()

    df = df[df["avg_depth"]<20].sort_values(by="avg_depth")
    print(len(df.index))
    fig = plt.figure(figsize=(20,10), dpi=300)
    plt.tight_layout()
    g = sns.barplot(data=df, x="MYID", y="avg_depth", color="blue")
    plt.xticks(rotation=90)
    plt.title("Low Coverage Mutations")
    plt.tight_layout()
    plt.savefig("low-coverage-mutations.pdf")
    plt.close()
    '''

                          
    # filtering
    genotypes = genotypes[genotypes['Waltz_total_t_depth'] >= min_depth]
    genotypes["AF"] = genotypes["Waltz_total_t_alt_count"] / genotypes["Waltz_total_t_depth"]


    # separate tumor and normal cells
    tumorCellsList = genotypes[(genotypes["MYID"].isin(truncalDriverMutations)) & (genotypes["AF"] >= min_AF_presence)]["Tumor_Sample_Barcode"].unique().tolist()
    tumorCellsCount = len(tumorCellsList)

    t = genotypes[genotypes["MYID"].isin(truncalDriverMutations)].copy()
    t["max_AF"] = t.groupby("Tumor_Sample_Barcode")["AF"].transform("max")
    normalCellsList = t[t["max_AF"]<=max_AF_absence]["Tumor_Sample_Barcode"].unique().tolist()
    normalCellsCount = len(normalCellsList)
    
    print(f"Total cells: {len(genotypes['Tumor_Sample_Barcode'].unique())}")
    print(f"Found {tumorCellsCount} cells with truncal driver mutations")
    print(f"Found {normalCellsCount} normal cells")

    tumorCells = genotypes[genotypes["Tumor_Sample_Barcode"].isin(tumorCellsList)].copy()
    normalCells = genotypes[genotypes["Tumor_Sample_Barcode"].isin(normalCellsList)].copy()

    #tumorCells.to_csv("tumor-cells.tsv", sep="\t", index=False)
    #normalCells.to_csv("normal-cells.tsv", sep="\t", index=False)


    # delete the original genotypes dataframe
    del genotypes


    '''
    # plot per-cell AF for truncal driver mutations
    df = tumorCells[tumorCells["MYID"].isin(truncalDriverMutations)].copy()
    sns.set_style("whitegrid")
    fig = plt.figure(dpi=300, figsize=(8,8))
    plt.tight_layout()
    g = sns.histplot(data=df, x="AF", bins=20, hue="MYID", multiple="dodge")
    sns.move_legend(g, "center")
    plt.title("Distribution of in-cell AFs for truncal driver mutations")
    plt.tight_layout()
    plt.savefig("truncal-driver-per-cell-AF-tumor.pdf")
    plt.close()

    df = normalCells[normalCells["MYID"].isin(truncalDriverMutations)].copy()
    sns.set_style("whitegrid")
    fig = plt.figure(dpi=300, figsize=(8,8))
    plt.tight_layout()
    g = sns.histplot(data=df, x="AF", bins=20, hue="MYID", multiple="dodge")
    sns.move_legend(g, "center")
    plt.title("Distribution of in-cell AFs for truncal driver mutations")
    plt.tight_layout()
    plt.savefig("truncal-driver-per-cell-AF-normal.pdf")
    plt.close()
    '''

    
    # average depth for each mutation within tumor cells and normal cells
    tumorCells["avg_depth_tumor"] = tumorCells.groupby("MYID")["Waltz_total_t_depth"].transform("mean")
    normalCells["avg_depth_normal"] = normalCells.groupby("MYID")["Waltz_total_t_depth"].transform("mean")

    '''
    # amplifications based on tumor and normal coverage
    t = tumorCells.drop_duplicates(subset=["MYID"])[["MYID", "avg_depth_tumor"]].copy()
    n = normalCells.drop_duplicates(subset=["MYID"])[["MYID", "avg_depth_normal"]].copy()

    t = pd.merge(t, n, on="MYID", how="left")
    t= t[~t["avg_depth_normal"].isna()]
    t["fc"] = t["avg_depth_tumor"] / t["avg_depth_normal"]
    t.to_csv("coverage-fc.tsv", sep="\t", index=False)


    sns.set_style("white")
    fig = plt.figure(dpi=300)
    plt.tight_layout()
    g = sns.scatterplot(data=t, x="MYID", y="fc", alpha=0.5)
    g.set(xticklabels=[])
    plt.xticks(rotation=90)
    plt.title("Fold change between average tumor and normal depths for each mutation")
    plt.tight_layout()
    plt.savefig("coverage-fc.pdf")
    plt.close()
    '''

    
    
    '''
    # allelic imbalance for FGFR1
    fgfr1snps = ["FGFR1_snp::8:38003787:T:C", "FGFR1_snp::8:38369982:C:G"]

    x = tumorCells[tumorCells["MYID"].isin(fgfr1snps)].copy()
    x["avg_af_tumor"] = x.groupby("MYID")["AF"].transform("mean")
    x = x.drop_duplicates(subset=["MYID"])[["MYID", "avg_af_tumor"]]
    print(x)

    y = normalCells[normalCells["MYID"].isin(fgfr1snps)].copy()
    y["avg_af_normal"] = y.groupby("MYID")["AF"].transform("mean")
    y = y.drop_duplicates(subset=["MYID"])[["MYID", "avg_af_normal"]]
    print(y)
    '''

    '''
    # distribution of per cell coverage fc for fgfr1
    # FGFR1
    snp = "cfDNA::8:38291950:G:A"
    # ELF3
    #snp = "ELF3_snp::1:201952574:T:C"
    #snp = "ELF3_snp::1:202092046:C:T"
    # SOX17
    #snp = "SOX17_snp::8:55049255:G:T"
    #snp = "SOX17_snp::8:55632762:C:T"
    # AKT3
    #snp = "AKT3_snp::1:243419429:T:A"
    #snp = "AKT3_snp::1:244601152:C:T"
    #ESR1
    #snp = "ESR1_snp::6:152446487:C:T"
    #snp = "ESR1_snp::6:152129077:T:C"

    t = tumorCells[tumorCells["MYID"] == snp].copy()

    n = normalCells[normalCells["MYID"] == snp].copy()
    n = n.drop_duplicates(subset=["MYID"])[["MYID", "avg_depth_normal"]]

    t = pd.merge(t, n, on="MYID", how="left")
    t["fc"] = t["Waltz_total_t_depth"] / t["avg_depth_normal"]
    
    sns.set_style("whitegrid")
    fig = plt.figure(dpi=300)
    plt.tight_layout()
    g = sns.histplot(data=t, x="fc")
    #g.set_xlim(0, 30)
    #sns.move_legend(g, "center")
    plt.title("Per cell fold change in depth for FGFR1")
    plt.tight_layout()
    plt.savefig("coverage-fc-cell-distribution.pdf")
    plt.close()
    '''
    

    
    # compute the cell fraction for each mutation, in tumor and in normal

    print("Computing cell fractions...")

    #genotypes["Present"] = ((genotypes["AF"] >= min_AF_presence) | ((genotypes["AF"]>=0.05) & (genotypes["Waltz_total_t_alt_count"]>=10)))

    tumorCells["Present"] = (tumorCells["AF"] >= min_AF_presence)                        
    tumorCells["PresentIn"] = tumorCells.groupby("MYID")["Present"].transform("sum")
    tumorCells["TumorCellFraction"] = tumorCells["PresentIn"] / tumorCellsCount

    normalCells["Present"] = (normalCells["AF"] >= min_AF_presence)                        
    normalCells["PresentIn"] = normalCells.groupby("MYID")["Present"].transform("sum")
    normalCells["NormalCellFraction"] = normalCells["PresentIn"] / normalCellsCount

    #tumorCells.to_csv("tumor-cells.tsv", sep="\t", index=False)
    #normalCells.to_csv("normal-cells.tsv", sep="\t", index=False)

    '''
    # distribution of af in tumor and normal cells
    sns.set_style("whitegrid")
    fig = plt.figure(dpi=300)
    plt.tight_layout()
    g = sns.ecdfplot(data=tumorCells[tumorCells["AF"]>0], x="AF")
    g.set_xlim(0, 0.2)
    #sns.move_legend(g, "center")
    plt.title("Distribution of AFs in tumor cells")
    plt.tight_layout()
    plt.savefig("af-distribution-tumor.pdf")
    plt.close()

    sns.set_style("whitegrid")
    fig = plt.figure(dpi=300)
    plt.tight_layout()
    g = sns.ecdfplot(data=normalCells[normalCells["AF"]>0], x="AF")
    g.set_xlim(0, 0.2)
    #sns.move_legend(g, "center")
    plt.title("Distribution of AFs in normal cells")
    plt.tight_layout()
    plt.savefig("af-distribution-normal.pdf")
    plt.close()
    '''

    '''
    # make a cf df
    cf = tumorCells.drop_duplicates(subset=["MYID"]).copy()
    cf = cf[["MYID", "avg_depth_tumor", "TumorCellFraction"]]

    x = normalCells.drop_duplicates(subset=["MYID"]).copy()
    x = x[["MYID", "avg_depth_normal", "NormalCellFraction"]]

    cf = pd.merge(cf, x, on="MYID", how="left")
    cf.to_csv("cf.tsv", sep="\t", index=False)

    # make a hotspots df
    hotspots = cf[(cf["MYID"].str.startswith("hotspot:") | cf["MYID"].str.startswith("ERBB2:") | cf["MYID"].str.startswith("PIK3CA:") | cf["MYID"].str.startswith("ESR1:") | cf["MYID"].str.startswith("TP53:") | cf["MYID"].str.startswith("ERBB3:"))]

    hotspots.to_csv("hotspots.tsv", sep="\t", index=False)
    '''

    '''
    # plot per cell af for mutations present in tissue but absent in single cells
    
    #m = "hotspot:ERBB2.p.L755S;Oncogenic:3A;HOTSPOT;17:37880220:T:C;;;Cancer"
    #m = "hotspot:TP53.p.D208N;LikelyOncogenic;HOTSPOT;17:7578227:C:T;;;Cancer"
    m = "hotspot:ESR1-IMPACT.p.Y537S;Oncogenic:3A;HOTSPOT;6:152419923:A:C;;;Cancer"
    
    df = normalCells[normalCells["MYID"] == m].copy()

    df.to_csv("t.tsv", sep="\t", index=False)

    sns.set_style("whitegrid")
    fig = plt.figure(dpi=300)
    plt.tight_layout()
    #g = sns.ecdfplot(data=df[df["AF"]>0], x="AF")
    g = sns.ecdfplot(data=df, x="AF")
    plt.title("Per cell AF distribution in normal cells for ESR1.Y537S")
    plt.tight_layout()
    plt.savefig("af-distribution-ESR1-normal.pdf")
    plt.close()
    '''

    
    # cluster hotspot mutations

    df = tumorCells[tumorCells["MYID"].str.startswith("hotspot:")].copy()
    df = df[df["TumorCellFraction"] >= 0.005]
    
    # remove BRCA2 mutation since it looks benign
    df = df[df["MYID"]!="hotspot:BRCA2-HRD.p.I423N;;;13:32906883:T:A;;;Cancer"]

    # make sure all mutations have sufficient coverage in all cells
    numMutations = df["MYID"].nunique()
    df["genotyped_muts"] = df.groupby("Tumor_Sample_Barcode")["MYID"].transform("nunique")
    df = df[df["genotyped_muts"] == numMutations]

    df = df[["Tumor_Sample_Barcode", "MYID", "Present", "PresentIn"]]
    df["PresentMutsNum"] = df.groupby("Tumor_Sample_Barcode")["Present"].transform("sum")
    df["I"] = df["Present"].astype(str)
    df["PresentMuts"] = df.groupby("Tumor_Sample_Barcode")["I"].transform(lambda x: ','.join(x))
    df = df.sort_values(by=["PresentMutsNum", "PresentMuts", "Tumor_Sample_Barcode", "PresentIn", "MYID"], ascending=[False, False, True, False, True])


    # check if ERBB2 mutations are in the same cells
    erbb21 = "hotspot:ERBB2.p.L755S;Oncogenic:3A;HOTSPOT;17:37880220:T:C;;;Cancer"
    erbb22 = "hotspot:ERBB2.p.D769Y;Oncogenic:3A;HOTSPOT;17:37880261:G:T;;;Cancer"
    t1 = df[((df["MYID"]==erbb21) & (df["Present"]==True))]["Tumor_Sample_Barcode"]
    t1 = set(t1)
    print(t1)
    t2 = df[((df["MYID"]==erbb22) & (df["Present"]==True))]["Tumor_Sample_Barcode"]
    t2 = set(t2)
    print(t2)
    print("intersection:", t1.intersection(t2))
    sys.exit(0)

    cols = df["Tumor_Sample_Barcode"].unique().tolist()
    rows = df["MYID"].unique().tolist()
    df = df.pivot(index="MYID", columns="Tumor_Sample_Barcode", values="Present").reindex(rows)[cols]
    
    df.to_csv("t.tsv", sep="\t", index=True)
    print(df.shape)

    sns.set_style("white")
    fig = plt.figure(dpi=300)
    plt.tight_layout()

    # distance mteric: hamming, rogerstanimoto, sokalmichener, cityblock - all work well, cosine, canberra
    r = sns.clustermap(df, method="complete", metric="hamming", col_cluster=False, row_cluster=False, xticklabels=False, yticklabels=True, cmap="Reds", figsize=(25, 10))
    r.ax_row_dendrogram.set_visible(False)
    r.ax_col_dendrogram.set_visible(False)
    r.cax.set_visible(False)
    r.ax_heatmap.set_xlabel("Cells", size=20)
    r.ax_heatmap.yaxis.set_ticks_position("left")

    # change label font size
    for label in r.ax_heatmap.get_yticklabels():
        label.set_size(20)

    plt.title("Hotspot/OncoKB mutations in individual cells")
    plt.tight_layout()
    plt.savefig("hotspot-clustering.pdf")
    plt.close()

    



    




    

   
    
    



    















    '''
    # get bulk CCFs
    bulkCCFs = pd.read_csv(bulkCCFsFile, sep="\t", index_col=False, keep_default_na=False, low_memory=False)

    bulkCCFs = bulkCCFs[bulkCCFs["Tumor_Sample_Barcode"] == sampleName]

    ccfs = pd.melt(bulkCCFs, id_vars="Tumor_Sample_Barcode", var_name='MYID', value_name='CCF_bulk')
    ccfs = ccfs.merge(ccfs_sc, on="MYID", how="left")
    ccfs = ccfs[~ccfs["CCF_sc"].isna()]
    ccfs["CCF_sc"] = ccfs["CCF_sc"].astype(float)
    ccfs["CCF_bulk"] = ccfs["CCF_bulk"].astype(float)

    ccfs.to_csv("ccfs.tsv", sep="\t", index=False)

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
    '''

        



    






    




if __name__ == '__main__':
    main()


