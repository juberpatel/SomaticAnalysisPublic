'''
Created on Jun 10, 2021

@author: patelj1
'''


'''
compute and plotAverages to answer
these questions:

1. What genes are mutated?
2. How many mutations per sample?
3. Average AF?
4. How many samples a gene is mutated in?


We will modify and improve this as we go on

'''

import os 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def makeAggDF(rowset, timepoint):
    # no of samples, no of mutations, avg af, handling rearrangements properly
    b = rowset.groupby("Treatment Arm").agg(
            NumMutations = pd.NamedAgg(column="AF", aggfunc=lambda s: len(s)),
            NumSamples = pd.NamedAgg(column="Sample", aggfunc=lambda s: len(np.unique(s))),
            NumSmallMutations = pd.NamedAgg(column="Hugo_Symbol", 
                                    aggfunc=lambda s: len([gene for gene in s if "__" not in gene])),
            SumAF = pd.NamedAgg(column="AF", aggfunc=lambda s: np.sum(s)),
            
            )
    
    b["Timepoint"] = timepoint
    b["AverageNumMutations"] = b["NumMutations"]/b["NumSamples"]
    b["AverageAF"] = b["SumAF"]/b["NumSmallMutations"]
    
    return b
    
def plotProgression(df, category):
    sns.set_theme(style="whitegrid")
    plt.figure(dpi=300) 
    plt.tight_layout()
    g = sns.catplot(data=df, x="Treatment Arm", hue="Progression", col=category,
        kind="count",
        height=6)
    
    g.despine(left=True)
    #g.set_ylabels(ylabel)
    
    
    
    plt.tight_layout()
    g.savefig(category + ".pdf")
        
    

def main():
    
    d = "/Users/patelj1/current/PromiseStudyJillian/analysis"
    mutationsFile = d + "/combined.tsv"
    
    os.chdir(d + "/patient-level")
    
    
    mutations = pd.read_csv(mutationsFile, sep="\t", index_col=False, keep_default_na=False, low_memory=False)
    
    # remove rearrangements
    # mutations = mutations[~mutations["Hugo_Symbol"].str.contains("__")]  
    # print(mutations["Hugo_Symbol"].value_counts())
    
    # get the sample lists
    #samples2 = patients.loc[(patients["Disease"] == "Breast") & (patients["CMO-ID"]!= ""), "CMO-ID"].tolist()
    
    pd.options.display.max_columns = None
    pd.options.display.max_rows = None
    
    patients = mutations.drop_duplicates(subset = ["Patient"])
    
    
    # progression status SBRT vs No SBRT based on
    # overall
    # Lung
    # Breast
    # Lung negative
    # Lung positive
    # Breast tnbc
    # Breat other
    
    
    df = patients.copy()
    plotProgression(df, "Disease")
    
    df = patients.copy().loc[(patients["Disease"]=="Lung") & (patients["Marker"]!="UNK")]
    df["Lung Driver"] = df["Marker"].apply(lambda x: "No" if x=="Negative" else "Yes")
    plotProgression(df, "Lung Driver")
    
    
    df = patients.copy().loc[(patients["Disease"]=="Breast") & (patients["Marker"]!="UNK")]
    df["Breast TNBC"] = "No"
    df.loc[((df["ER"]=="-") & (df["PR"]=="-") & (df["HER2"]=="-")), "Breast TNBC"] = "Yes"
    plotProgression(df, "Breast TNBC")
    
    
    
    '''
    b = makeAggDF(baseline, "Baseline")
    e = makeAggDF(eightWeeks, "8-week")
    p = makeAggDF(progression, "Progression")
    
    df = b.append(e).append(p)
    df = df.reset_index()
    #print(df)
    
    # no of samples, no of mutations, avg af, handling rearrangements properly
    plotAverages(df, "NumSamples", "Number of Samples")
    plotAverages(df, "AverageNumMutations", "Avg. Mutations per Sample")
    plotAverages(df, "AverageAF", "Average AF")
    '''
    
    '''
    # add slight AF for rearrangements to make them count
    baseline.loc[baseline["Hugo_Symbol"].str.contains("__"), "AF"] = 0.0001
    eightWeeks.loc[eightWeeks["Hugo_Symbol"].str.contains("__"), "AF"] = 0.0001
    progression.loc[progression["Hugo_Symbol"].str.contains("__"), "AF"] = 0.0001
    '''
    '''
    # plot gene counts
    b = makeGeneCountDF(baseline, "Baseline")
    e = makeGeneCountDF(eightWeeks, "8-week")
    p = makeGeneCountDF(progression, "Progression")
    
    
    
    df = b.append(e).append(p)
    df = df.reset_index()
    
    print(df)
    
    plotGeneCounts(df, "Frequency", "Frequency")
    '''
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    '''
    # get the mutated genes
    genes1 = mutations.loc[(mutations["Patient"].isin(samples1)) 
                           & (mutations["SampleType"]=="cfDNA") ,"Hugo_Symbol"
                           ].value_counts().rename_axis('Gene').reset_index(name='# Samples')
    '''
    
    
    
    
    
   # print(genes1)

    
    
    
    
    
    
    
    
















if __name__ == '__main__':
    main()
    
    
    
    
    