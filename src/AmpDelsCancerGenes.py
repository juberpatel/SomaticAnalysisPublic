'''
Created on Nov 3, 2022

@author: Juber Patel

'''

import pandas as pd
import os
import seaborn as sb
import copy
import matplotlib.pyplot as plt
import statistics
from collections import OrderedDict
import sys



def getLCN(sample, row, copyNumbers):
    
    x = copyNumbers[(copyNumbers["ID"]==sample) & (copyNumbers["chrom"]==row["Chromosome"]) 
                    & (copyNumbers["loc.start"]<=row["Start"]) & (copyNumbers["loc.end"]>=row["Start"])]
    
    
    x = x.iloc[:1]["tcn.em"]
    
    tcn = 2
    
    if len(x) != 0:
        tcn = x.values[0] 
    
    return tcn


def getInterestingGenes(row, geneCoordinates, cancerGenes):
    
    x = geneCoordinates[(row["chrom"]==geneCoordinates["Chromosome"]) & (row["loc.start"]<=geneCoordinates["Start"])
                               & (row["loc.end"]>=geneCoordinates["Start"])]
    
    x = x.drop_duplicates(subset=["Name"])
    x = x[x["Name"].isin(cancerGenes)]
    
    return x
    
    #genes = x["Name"].asList()
    
    '''
    for gene in genes:
        if gene in cancerGenes:
            ampDelGenes.append({"Gene": gene, "Chromosome": })
            
        elif gene + "-HRD" in cancerGenes:
            ampDelGenes.add(gene + "-HRD")
            
        elif gene + "-IMPACT" in cancerGenes:
            ampDelGenes.add(gene + "-IMPACT")
    '''       
        
    
    
    


def main():
    
    
    patient = "MSK-AB-0002"
    d = "/Users/patelj1/current/Autopsy/analysis/filtered-new/" + patient
    os.chdir(d)

    patientName = patient
    badSamplesFile = patient + "-bad-samples.txt"
    copyNumberFile = patient + "-copy-numbers.txt"
    
    
    geneCoordinatesFile = "/Users/patelj1/resources/gene-list/juber-hg19-gene-list.bed"
    geneListFile = "/Users/patelj1/current/workspace/SomaticAnalysis/cancer-HRD-impact505-genes.txt"
    sampleLabelsFile = "/Users/patelj1/current/Autopsy/analysis/labels.txt"
    
    
    # load cancer genes
    cancerGenes = set()
    with open(geneListFile) as f:
        for line in f:
            gene = line.split()[0]
            cancerGenes.add(gene)
            
   
    # load gene coordinates
    geneCoordinates = pd.read_csv(geneCoordinatesFile, sep="\t", index_col=False, keep_default_na=False, low_memory=False)
    geneCoordinates["Chromosome"] = geneCoordinates["Chromosome"].astype(str)
    #geneCoordinates["Start"] = geneCoordinates["Start"].astype(int)
    geneCoordinates["Name"] = geneCoordinates["Name"].apply(lambda x: x.split("_")[0])
    geneCoordinates["Name"] = geneCoordinates["Name"].apply(lambda x: x + "-HRD" if (x + "-HRD") in cancerGenes else x)
    geneCoordinates["Name"] = geneCoordinates["Name"].apply(lambda x: x + "-IMPACT" if (x + "-IMPACT") in cancerGenes else x)
    
    
    amplifiedGenes = pd.DataFrame()
    deletedGenes = pd.DataFrame()
    
    # go through copy numbers file and find out the segments that have amplification or homdel
    copyNumbers = pd.read_csv(copyNumberFile, sep="\t", index_col=False, keep_default_na=False, low_memory=False)
    copyNumbers["ID"] = copyNumbers["ID"].apply(lambda x: x.split("_")[0])
    copyNumbers["chrom"] = copyNumbers["chrom"].astype(str)
    #copyNumbers["loc.start"] = copyNumbers["loc.start"].astype(int)
    copyNumbers.loc[copyNumbers["chrom"]=="23", "chrom"] = "X"
    
    # drop the specified bad samples  
    badSamples = []
    with open(badSamplesFile) as f:
        for line in f:
            badSamples.append(line.split()[0])
    
    print("Removing Samples:")
    print(badSamples)
    copyNumbers = copyNumbers[~(copyNumbers["ID"].isin(badSamples))]
    goodSamples = copyNumbers["ID"].unique().tolist()
    
    
    
    amplifiedSegments = copyNumbers[copyNumbers["tcn.em"]>=7]
    deletedSegments = copyNumbers[copyNumbers["tcn.em"]==0]
    
    for index, row in amplifiedSegments.iterrows():
        x = getInterestingGenes(row, geneCoordinates, cancerGenes)
        amplifiedGenes = pd.concat([amplifiedGenes, x])

    for index, row in deletedSegments.iterrows():
        x = getInterestingGenes(row, geneCoordinates, cancerGenes)
        deletedGenes = pd.concat([deletedGenes, x])
    
    amplifiedGenes = amplifiedGenes.drop_duplicates(subset=["Name"])
    deletedGenes = deletedGenes.drop_duplicates(subset=["Name"])
    
    print(len(amplifiedGenes))
    
    if len(amplifiedGenes) > 0:
        print(sorted(amplifiedGenes["Name"].tolist()))
    
    print(len(deletedGenes))
    
    if len(deletedGenes) > 0:
        print(sorted(deletedGenes["Name"].tolist()))
    
    
    ampDelGenes = pd.concat([amplifiedGenes, deletedGenes])
    ampDelGenes = ampDelGenes.drop_duplicates(subset=["Name"])
    
    if len(ampDelGenes) == 0:
        sys.exit()
    
    ampDelGenes["Name"] = ampDelGenes["Name"] + ":" + ampDelGenes["Chromosome"] + ":" + ampDelGenes["Start"].astype(str)
    
    df = pd.DataFrame(columns=["Sample", "Gene", "tcn"])
    
    # find out tcn for each of the interesting genes in all samples
    for sample in goodSamples:
        for index, row in ampDelGenes.iterrows():
            tcn = getLCN(sample, row, copyNumbers)
            
            df.loc[len(df.index)] = [sample, row["Name"], tcn]
            
    
    # use provided sample labels
    sampleLabels = {}
    with open(sampleLabelsFile) as f:
        for line in f:
            tokens = line.split()
            sampleLabels[tokens[0]] = tokens[1]
    
    df = df.replace({"Sample": sampleLabels})
    
    
    # make the CCF matrix 
    ampDelMatrix = df.pivot_table(index="Gene", columns="Sample", values="tcn", fill_value=2)
    
    
    # sort on number of samples a mutation is present in, then the exact samples it is present in,
    # then the sum of CCF values in those samples
    
    amplifiedIn = df.copy()
    amplifiedIn.loc[amplifiedIn["tcn"]<7, "tcn"] = 0
    amplifiedIn = amplifiedIn.pivot_table(index="Gene", columns="Sample", values="tcn", fill_value=0)
    amplifiedIn = amplifiedIn.astype(bool).replace({True:"1", False:"0"}).apply("".join, axis=1)
    
    
    
    deletedIn = df.copy()
    deletedIn.loc[deletedIn["tcn"]>0, "tcn"] = 0
    deletedIn = deletedIn.pivot_table(index="Gene", columns="Sample", values="tcn", fill_value=0)
    deletedIn = deletedIn.astype(bool).replace({True:"1", False:"0"}).apply("".join, axis=1)
     
    sums = ampDelMatrix.sum(axis=1)

    # add these columns AFTER calculating them outside of the dataframe
    ampDelMatrix["sums"] = sums
    ampDelMatrix["amplifiedIn"] = amplifiedIn
    ampDelMatrix["deletedIn"] = deletedIn
    
    
    # do the sorting
    #ampDelMatrix.sort_values(by=["rank", "sums"], ascending=False, inplace=True)
    ampDelMatrix.sort_values(by=["sums"], ascending=False, inplace=True)
    
    
    '''
    # get cluster number for each mutation
    count = 1
    mutationClusters = {}
    groups = ampDelMatrix.groupby(["presentIn"], sort=False).groups
    for group, labels in groups.items():
        print(group + "\t" + str(len(labels)))
        for label in labels:
            mutationClusters[label] = count
            
        count += 1
    
    # add cluster number to df
    df["Cluster"] = df.apply(lambda row: mutationClusters.get(row["MYID"], -1), axis=1)
    '''
    
    # save matrix and mutations to files
    ampDelMatrix = ampDelMatrix.drop(["sums", "amplifiedIn", "deletedIn"], axis=1)
    t = ampDelMatrix.T
    ampDelMatrixFile = patientName + "-AmpDelMatrix.tsv"
    t.to_csv(ampDelMatrixFile, sep="\t")
    
    '''
    clusteredFile = mutationsFile.replace(".tsv", "-clustered.tsv")
    df.to_csv(clusteredFile, sep="\t", index=False)
    '''
    
    
    
    print("Plotting AmpDel heatmap...")
    AmpDelMatrix = pd.read_csv(ampDelMatrixFile, sep="\t", index_col=0, low_memory=False)
   
    '''
    # get cluster numbers
    clusterNumbers = []
    for mut in AmpDelMatrix.columns:
        clusterNumbers.append(mutationClusters[mut])
        
    print(clusterNumbers)
    '''
    
    fig = plt.figure(dpi=300)
    plt.tight_layout()
    
    
    c = sb.color_palette("Reds", as_cmap=True)
    c = copy.copy(c)
    c.set_under("white")
    
    '''
    # make cluster marker colors
    #columnPalette = sb.color_palette(["black", "white"])
    columnPalette = sb.color_palette()
    columnColors = {}
    for mut in mutationClusters:
        indicator = mutationClusters[mut] % 2
        columnColors[mut] = columnPalette[indicator]
        
    columnColors = pd.Series(columnColors)
    '''
    
    r = sb.clustermap(AmpDelMatrix, method="complete", metric="canberra", col_cluster=False, 
                      cmap=c, vmin=1, vmax=10, cbar_kws={'label': 'tcn'}, figsize=(20, 10),
                      xticklabels=True, yticklabels=True)
    r.ax_row_dendrogram.set_visible(False)
    
    '''
    labels = []
    for tick in r.ax_heatmap.xaxis.get_major_ticks():
        label = tick.label1
        label.set_size(15)
        text = label.get_text()
        name = text.split(";")[0]
        #name = text
        
        if "HOTSPOT" in text:   
            label.set_visible(True)
            label.set_color("red")
        elif "Oncogenic" in text:
            label.set_visible(True)
            label.set_color("blue")
        else:
            label.set_visible(False)
            tick.tick1line.set_visible(False)
            tick.label1.set_visible(False)
            
            
        labels.append(name)
        
    
    r.ax_heatmap.set_xticklabels(labels)
    '''
    
    
    
    plt.savefig(ampDelMatrixFile.replace(".tsv", ".pdf"))
    
    
    
    
    
    
    
    
    
    





if __name__ == '__main__':
    main()
