'''
Created on November 30, 2023

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
    #bulkCCFsFile = "MSK-AB-0002-absolute-filtered-annotated-ccfs.txt"
    bulkCCFsFile = "MSK-AB-0002-ABSOLUTE-oncokb-clustered.tsv"
    bulkGenotypingFile = "../../../MSK-AB-0002.genotyped_mutations.jlab_and_cmo.filtered.maf-exac-annoated.tsv"
    readCountsFile = "read-counts.txt"


    truncalDriverMutations = ["hotspot:CDH1.p.W156Tfs*60;LikelyOncogenic;;16:68842402:G:GAC;;;Cancer", "hotspot:TP53.p.N239D;LikelyOncogenic;HOTSPOT;17:7577566:T:C;;;Cancer"]

    sampleName = "T23_1-left-top-back-skin-met-1"

    min_depth = 15
    min_AF_presence = 0.025
    max_AF_absence = 0.02

    genotypes = pd.read_csv(genotypingFile, sep="\t", index_col=False, keep_default_na=False, low_memory=False)


    bulkCCFs = pd.read_csv(bulkCCFsFile, sep="\t", index_col=False, keep_default_na=False, low_memory=False)
    bulkCCFs = bulkCCFs[["Tumor_Sample_Barcode", "ID", "MYID", "Cluster", "Cancer_Cell_Fraction"]]

    bulkGenotypes = pd.read_csv(bulkGenotypingFile, sep="\t", index_col=False, keep_default_na=False, low_memory=False)
    bulkGenotypes = bulkGenotypes[["Tumor_Sample_Barcode", "ID", "Hugo_Symbol", "HGVSp_Short", "Variant_Classification", "Variant_Type", "TUMOR_MAF", "TUMOR_DP"]]
    bulkGenotypes = bulkGenotypes[bulkGenotypes["Tumor_Sample_Barcode"] == "MSK-AB-0002-T23-1"]


    readCounts = pd.read_csv(readCountsFile, sep="\t", index_col=False, keep_default_na=False, low_memory=False)
    readCounts = readCounts[["ID",  "TotalOnTarget"]]
    readCounts["ID"] = readCounts["ID"].str.split(".bam").str[0]
    readCounts["ID"] = readCounts["ID"].str.replace("_", "-")
    readCounts = readCounts.rename(columns={"ID": "Tumor_Sample_Barcode", "TotalOnTarget": "CellOnTargetReads"})
    readCounts["TotalReads"] = readCounts["CellOnTargetReads"].sum()
    readCounts["CellOnTargetReadsFraction"] = readCounts["CellOnTargetReads"] / readCounts["TotalReads"]
    

    os.chdir("plots")

    sns.set_style("whitegrid")

    # add read counts
    genotypes = pd.merge(genotypes, readCounts, on="Tumor_Sample_Barcode", how="left")

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
    # plot distribution of reads among tumor and normal cells
    tumorCells["CellType"] = "Tumor"
    normalCells["CellType"] = "Normal"
    df = pd.concat([tumorCells, normalCells])
    df = df[['Tumor_Sample_Barcode', 'CellOnTargetReads', 'CellType']]
    df = df.drop_duplicates(subset=["Tumor_Sample_Barcode"])

    fig = plt.figure(dpi=300)
    plt.tight_layout()

    g = sns.boxplot(data=df, x="CellType", y="CellOnTargetReads", showfliers=False, color="white")
    sns.stripplot(data=df, x="CellType", y="CellOnTargetReads", dodge=True, ax=g, alpha=0.7, color="red", jitter=0.4, size=3)
    plt.ticklabel_format(style='plain', axis='y')

    plt.title("Reads per cell")
    plt.tight_layout()
    plt.savefig("reads-per-cell.pdf")
    plt.close()
    '''


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

    '''
    # average depth for each mutation within tumor cells
    tumorCells["avg_depth_tumor"] = tumorCells.groupby("MYID")["Waltz_total_t_depth"].transform("mean")
    normalCells["avg_depth_normal"] = normalCells.groupby("MYID")["Waltz_total_t_depth"].transform("mean")

    
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


    '''
    # make a heatmap of fold changes for interesting snps

    snps = ["cfDNA::8:38291950:G:A", "AKT3_snp::1:243419429:T:A"]

    df = tumorCells[tumorCells["MYID"].isin(snps)].copy()

    n = normalCells[normalCells["MYID"].isin(snps)].copy()
    n = n.drop_duplicates(subset=["MYID"])[["MYID", "avg_depth_normal"]]

    df = pd.merge(df, n, on="MYID", how="left")
    df["fc"] = df["Waltz_total_t_depth"] / df["avg_depth_normal"]
    df.loc[df["fc"]>5, "fc"] = 5
    # change the snp name slightly
    df.loc[df["MYID"]=="cfDNA::8:38291950:G:A", "MYID"] = "FGFR1_cfDNA::8:38291950:G:A"


    # make sure all mutations have sufficient coverage in all cells
    numMutations = df["MYID"].nunique()
    df["genotyped_muts"] = df.groupby("Tumor_Sample_Barcode")["MYID"].transform("nunique")
    df = df[df["genotyped_muts"] == numMutations]


    #df = df[["Tumor_Sample_Barcode", "MYID", "Present", "PresentIn"]]
    #df["PresentMutsNum"] = df.groupby("Tumor_Sample_Barcode")["Present"].transform("sum")
    #df["I"] = df["Present"].astype(str)
    #df["PresentMuts"] = df.groupby("Tumor_Sample_Barcode")["I"].transform(lambda x: ','.join(x))
    #df = df.sort_values(by=["PresentMutsNum", "PresentMuts", "Tumor_Sample_Barcode", "PresentIn", "MYID"], ascending=[False, False, True, False, True])

    df["avg_fc"] = df.groupby("Tumor_Sample_Barcode")["fc"].transform("mean")
    df = df.sort_values(by=["avg_fc", "Tumor_Sample_Barcode", "MYID"])



    cols = df["Tumor_Sample_Barcode"].unique().tolist()
    rows = df["MYID"].unique().tolist()
    df = df.pivot(index="MYID", columns="Tumor_Sample_Barcode", values="fc").reindex(rows)[cols]
    
    df.to_csv("t.tsv", sep="\t", index=True)
    print(df.shape)

    sns.set_style("white")
    fig = plt.figure(dpi=300)
    plt.tight_layout()

    # RdBu_r, vlag
    g = sns.heatmap(df, center=1, cmap="RdBu_r")
    g.set(xlabel="Cells")
    g.set(xticklabels=[])

    plt.title("Depth fc in individual cells")
    plt.tight_layout()
    plt.savefig("coverage-fc-heatmap.pdf")
    plt.close()
    '''
    

    
    # compute the cell fraction for each mutation, in tumor and in normal

    print("Computing cell fractions...")

    #genotypes["Present"] = ((genotypes["AF"] >= min_AF_presence) | ((genotypes["AF"]>=0.05) & (genotypes["Waltz_total_t_alt_count"]>=10)))

    tumorCells["Present"] = (tumorCells["AF"] >= min_AF_presence)                        
    tumorCells["PresentIn"] = tumorCells.groupby("MYID")["Present"].transform("sum")
    tumorCells["TumorCellFraction"] = tumorCells["PresentIn"] / tumorCellsCount

    tumorCells["CollectiveDepth"] = tumorCells.groupby("MYID")["Waltz_total_t_depth"].transform("sum")
    tumorCells["CollectiveAltCount"] = tumorCells.groupby("MYID")["Waltz_total_t_alt_count"].transform("sum")
    tumorCells["CollectiveAF"] = tumorCells["CollectiveAltCount"] / tumorCells["CollectiveDepth"]


    normalCells["Present"] = (normalCells["AF"] >= min_AF_presence)                        
    normalCells["PresentIn"] = normalCells.groupby("MYID")["Present"].transform("sum")
    normalCells["NormalCellFraction"] = normalCells["PresentIn"] / normalCellsCount

    normalCells["CollectiveDepth"] = normalCells.groupby("MYID")["Waltz_total_t_depth"].transform("sum")
    normalCells["CollectiveAltCount"] = normalCells.groupby("MYID")["Waltz_total_t_alt_count"].transform("sum")
    normalCells["CollectiveAF"] = normalCells["CollectiveAltCount"] / normalCells["CollectiveDepth"]


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
    #sys.exit(0)
    

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


    '''
    # plot TumorCellFraction, CollectiveAF, bulk AF and bulk depth for clusters not present in bulk T23-1

    bulkCCFs["Cluster"] = bulkCCFs["Cluster"].astype(int)
    allClusters = set(bulkCCFs["Cluster"].unique())
    presentClusters = set(bulkCCFs[(bulkCCFs["Tumor_Sample_Barcode"]==sampleName)]["Cluster"].unique())
    absentClusters = allClusters - presentClusters

    interestingMutations = bulkCCFs[bulkCCFs["Cluster"].isin(absentClusters)]["MYID"].unique().tolist()
    print(len(interestingMutations))
    
    
    x = bulkCCFs[bulkCCFs["MYID"].isin(interestingMutations)].copy()
    x = x.drop_duplicates(subset=["MYID"])[["MYID", "ID", "Cluster"]]
    x.to_csv("x.tsv", sep="\t", index=False)
    
    
    tumorCells["MYID"] = tumorCells["MYID"].str.split(":").str[1:].str.join(":")
    df = tumorCells[tumorCells["MYID"].isin(interestingMutations)].copy()
    df = df.drop_duplicates(subset=["MYID"])
    df = pd.merge(df, x, on="MYID", how="left")
    interestingMutations = df["MYID"].tolist()
    print(len(interestingMutations))
    #df.to_csv("t.tsv", sep="\t", index=False)
    
    # plot TumorCellFractions
    sns.set_style("white")
    fig = plt.figure(dpi=300)
    plt.tight_layout()

    #boxprops={'alpha': 0.2},
    g = sns.boxplot(data=df, x="Cluster", y="TumorCellFraction", showfliers=False, color="white")
    sns.stripplot(data=df, x="Cluster", y="TumorCellFraction", dodge=True, ax=g, alpha=0.7, color="red", jitter=0.4, size=3)
    #sns.swarmplot(data=df, x="Cluster", y="TumorCellFraction", ax=g, color="red", size=3, alpha=0.5)
    
    plt.title("SC TumorCellFraction for clusters not present in bulk T23-1")
    plt.tight_layout()
    plt.savefig("absent-clusters-tumor-fractions-0.025.pdf")
    plt.close()

    # plot CollectiveAFs
    sns.set_style("white")
    fig = plt.figure(dpi=300)
    plt.tight_layout()

    #boxprops={'alpha': 0.2},
    g = sns.boxplot(data=df, x="Cluster", y="CollectiveAF", showfliers=False, color="white")
    sns.stripplot(data=df, x="Cluster", y="CollectiveAF", dodge=True, ax=g, alpha=0.7, color="red", jitter=0.4, size=3)
    #sns.swarmplot(data=df, x="Cluster", y="CollectiveAF", ax=g, color="red", size=3)
    
    plt.title("SC CollectiveAF for clusters not present in bulk T23-1")
    plt.tight_layout()
    plt.savefig("absent-clusters-collective-af.pdf")
    plt.close()


    bulkGenotypes = pd.merge(bulkGenotypes, x, on="ID", how="left")
    df = bulkGenotypes[bulkGenotypes["MYID"].isin(interestingMutations)].copy()
    df["Cluster"] = df["Cluster"].astype(int)
    
    # plot bulk AFs
    sns.set_style("white")
    fig = plt.figure(dpi=300)
    plt.tight_layout()

    #boxprops={'alpha': 0.2},
    g = sns.boxplot(data=df, x="Cluster", y="TUMOR_MAF", showfliers=False, color="white")
    sns.stripplot(data=df, x="Cluster", y="TUMOR_MAF", dodge=True, ax=g, alpha=0.7, color="red", jitter=0.4, size=3)
    #sns.swarmplot(data=df, x="Cluster", y="TUMOR_MAF", ax=g, color="red", size=3)
    
    plt.title("Bulk AFs for clusters not present in bulk T23-1")
    plt.tight_layout()
    plt.savefig("absent-clusters-bulk-af.pdf")
    plt.close()

    # plot bulk depths
    sns.set_style("white")
    fig = plt.figure(dpi=300)
    plt.tight_layout()

    #boxprops={'alpha': 0.2},
    g = sns.boxplot(data=df, x="Cluster", y="TUMOR_DP", showfliers=False, color="white")
    sns.stripplot(data=df, x="Cluster", y="TUMOR_DP", dodge=True, ax=g, alpha=0.7, color="red", jitter=0.4, size=3)
    #sns.swarmplot(data=df, x="Cluster", y="TUMOR_DP", ax=g, color="red", size=3)
    
    plt.title("Bulk depths for clusters not present in bulk T23-1")
    plt.tight_layout()
    plt.savefig("absent-clusters-bulk-depth.pdf")
    plt.close()
    '''

    '''
    # plot sc tcf of interesting 'absent' clusters (11, 13, 20) vs. ccf in bulk T43-1
    bulkCCFs["Cluster"] = bulkCCFs["Cluster"].astype(int)
    x = bulkCCFs[bulkCCFs["Tumor_Sample_Barcode"]=="T43_1-right-lower-lobe-lung-LN"].copy()
    x = x[x["Cluster"].isin([11, 13, 20])]
    x.to_csv("x.tsv", sep="\t", index=False)
    
    # get single cell MYID in line with bulk MYID
    tumorCells["MYID"] = tumorCells["MYID"].str.split(":").str[1:].str.join(":")
    df = tumorCells[tumorCells["MYID"].isin(x["MYID"])].copy()
    df = df.drop_duplicates(subset=["MYID"])
    df = pd.merge(df, x, on="MYID", how="left")
    interestingMutations = df["MYID"].tolist()
    print(len(interestingMutations))
    df.to_csv("t.tsv", sep="\t", index=False)

    
    # plot CCFs vs. TumorCellFractions
    sns.set_style("whitegrid")
    fig = plt.figure(dpi=300)
    plt.tight_layout()

    g = sns.scatterplot(data=df, x="Cancer_Cell_Fraction", y="TumorCellFraction", hue="Cluster", palette="tab10")
    g.set(xlabel="T43-1 Bulk Cancer Cell Fraction", ylabel="T23-1 SC Tumor Cell Fraction")
    
    plt.title("T43-1 Bulk CCF vs. T23-1 SC TumorCellFraction")
    plt.tight_layout()
    plt.savefig("t43-1-bulk-ccf-t23-1-sc-tcf.pdf")
    plt.close()
    '''

    


    '''
    # plot bulk CCF vs. TumorCellFraction for truncal clonal mutations

    bulkCCFs["Cluster"] = bulkCCFs["Cluster"].astype(int)
    bulkCCFs["avg_CCF"] = bulkCCFs.groupby("MYID")["Cancer_Cell_Fraction"].transform("mean")
    x = bulkCCFs[bulkCCFs["Tumor_Sample_Barcode"]==sampleName].copy()
    x = x[((x["Cluster"]==1) & (x["avg_CCF"]>0.75))]
    x.to_csv("x.tsv", sep="\t", index=False)
    
    # get single cell MYID in line with bulk MYID
    tumorCells["MYID"] = tumorCells["MYID"].str.split(":").str[1:].str.join(":")
    df = tumorCells[tumorCells["MYID"].isin(x["MYID"])].copy()
    df = df.drop_duplicates(subset=["MYID"])
    df = pd.merge(df, x, on="MYID", how="left")
    interestingMutations = df["MYID"].tolist()
    print(len(interestingMutations))
    #df.to_csv("t.tsv", sep="\t", index=False)

    
    # plot CCFs vs. TumorCellFractions
    sns.set_style("whitegrid")
    fig = plt.figure(dpi=300)
    plt.tight_layout()

    g = sns.scatterplot(data=df, x="Cancer_Cell_Fraction", y="TumorCellFraction")
    g.set(xlabel="Bulk Cancer Cell Fraction", ylabel="SC Tumor Cell Fraction")
    
    plt.title("Bulk CCF vs. SC TumorCellFraction")
    plt.tight_layout()
    plt.savefig("bulk-ccf-sc-tcf-0.15.pdf")
    plt.close()
    '''


    '''
    # plot clusters covered in the panel and their mutation counts

    bulkCCFs["Cluster"] = bulkCCFs["Cluster"].astype(int)
    x = bulkCCFs.drop_duplicates(subset=["MYID"]).copy()
    x = x[x["Cluster"]!=-1]

    # get single cell MYID in line with bulk MYID
    tumorCells["MYID"] = tumorCells["MYID"].str.split(":").str[1:].str.join(":")
    df = tumorCells[tumorCells["MYID"].isin(x["MYID"])].copy()
    df = df.drop_duplicates(subset=["MYID"])
    df = pd.merge(df, x, on="MYID", how="left")
    df["ClusterCount"] = df.groupby("Cluster")["MYID"].transform("nunique")
    interestingMutations = df["MYID"].tolist()
    print(len(interestingMutations))
    df.to_csv("t.tsv", sep="\t", index=False)
    
    
    # plot the clusters covered in the panel
    sns.set_style("whitegrid")
    fig = plt.figure(dpi=300)
    plt.tight_layout()

    g = sns.barplot(data=df, x="Cluster", y="ClusterCount", color="blue", errorbar=None)
    g.set(ylabel="Mutation Count")
    
    plt.title("Clusters covered in the panel")
    plt.tight_layout()
    plt.savefig("covered-clusters.pdf")
    plt.close()
    '''


    # identify and plot clusters present in sc T23-1 and plot them against individual cells

    bulkCCFs["Cluster"] = bulkCCFs["Cluster"].astype(int)
    x = bulkCCFs.drop_duplicates(subset=["MYID"]).copy()
    x = x.drop("Tumor_Sample_Barcode", axis=1)
    x = x[x["Cluster"]!=-1]
    # remove clusters 7 because it has only 1 mutation which is not covred well
    x = x[x["Cluster"]!=7]
    #x = x[x["Cluster"]!=33]

    # get single cell MYID in line with bulk MYID
    tumorCells["MYID"] = tumorCells["MYID"].str.split(":").str[1:].str.join(":")
    df = tumorCells[tumorCells["MYID"].isin(x["MYID"])].copy()
    df = pd.merge(df, x, on="MYID", how="left")
    df["ClusterCount"] = df.groupby("Cluster")["MYID"].transform("nunique")
    df["ClusterThreshold"] = df["ClusterCount"] / 10.0
    df.loc[df["ClusterThreshold"]<1.0, "ClusterThreshold"] = 1.0

    df["ClusterPresentMutations"] = df.groupby(["Tumor_Sample_Barcode", "Cluster"])["Present"].transform("sum")
    df["ClusterPresent"] = (df["ClusterPresentMutations"] >= df["ClusterThreshold"])
    #df["ClusterPresent"] = True
    df = df[['Tumor_Sample_Barcode', 'Cluster', 'ClusterPresent', 'ClusterPresentMutations', 'ClusterThreshold', 'ClusterCount']]
    df = df.drop_duplicates(subset=["Tumor_Sample_Barcode", "Cluster"])
    
    # remove cells that do not have 32 clusters
    df["NumClusters"] = df.groupby("Tumor_Sample_Barcode")["Cluster"].transform("nunique")
    df = df[df["NumClusters"] == 32]
    df = df.drop("NumClusters", axis=1)

    # cluster the cells and clusters
    #df = df[["Tumor_Sample_Barcode", "MYID", "Present", "PresentIn"]]
    df["PresentClustersNum"] = df.groupby("Tumor_Sample_Barcode")["ClusterPresent"].transform("sum")
    df["I"] = df["ClusterPresent"].astype(str)
    df["PresentClusters"] = df.groupby("Tumor_Sample_Barcode")["I"].transform(lambda x: ','.join(x))
    df["ClusterPresentIn"] = df.groupby("Cluster")["ClusterPresent"].transform("sum")
    num_cells = df["Tumor_Sample_Barcode"].nunique()
    df["TumorCellFraction"] = df["ClusterPresentIn"] / num_cells
    df = df[df["TumorCellFraction"] >= 0.02]
    #df = df[df["Cluster"].isin([1, 2, 5, 8, 11, 20, 13, 21, 33, 31, 10, 9])]
    df = df.sort_values(by=["PresentClustersNum", "PresentClusters", "Tumor_Sample_Barcode", "ClusterPresentIn", "Cluster"], ascending=[False, False, True, False, True])
    df.to_csv("clusters-in-cells.tsv", sep="\t", index=False)


    cols = df["Tumor_Sample_Barcode"].unique().tolist()
    rows = df["Cluster"].unique().tolist()
    df = df.pivot(index="Cluster", columns="Tumor_Sample_Barcode", values="ClusterPresent").reindex(rows)[cols]
    #df.to_csv("t.tsv", sep="\t", index=True)
    print(df.shape)


    sns.set_style("white")
    fig = plt.figure(dpi=300, figsize=(8,6))
    plt.tight_layout()

    # RdBu_r, vlag
    g = sns.heatmap(df, cmap="Reds", cbar=False)
    g.set(xlabel="Cells")
    g.set(xticklabels=[])

    plt.yticks(rotation=0)
    plt.title("Mutation clusters in individual cells")
    plt.tight_layout()
    plt.savefig("clusters-in-cells-filtered.pdf")
    plt.close()

    

    
    
    



    

   
    
    



    















    






    




if __name__ == '__main__':
    main()


