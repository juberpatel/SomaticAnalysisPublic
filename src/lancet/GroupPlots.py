'''
Created on Jun 10, 2021

@author: Juber Patel

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
from datetime import date


def makeGeneCountDF(rowset, timepoint):
    
    b = rowset.groupby(["Treatment Arm", "Hugo_Symbol"]).agg(
            Frequency = pd.NamedAgg(column="AF", aggfunc=lambda s: np.count_nonzero(s)),
            )
    
    b = b.reset_index()
    b = b.groupby(["Treatment Arm"]).apply(lambda g: g.nlargest(10, "Frequency", keep="all"))
    
     
    b["Timepoint"] = timepoint
    b = b.reset_index(drop=True)
    
    return b


def plotGeneCounts(df, yaxis, ylabel):
    sns.set_theme(style="whitegrid")
    plt.figure(dpi=300) 
    plt.tight_layout()
    g = sns.catplot(data=df, x="Treatment Arm", y=yaxis, hue="Hugo_Symbol", col="Timepoint",
        kind="bar",
        height=6)
    
    g.despine(left=True)
    g.set_ylabels(ylabel)
    
    
    
    plt.tight_layout()
    g.savefig(ylabel + ".pdf")
    

#############################################################



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
    
def plotAverages(df, yaxis, ylabel):
    sns.set_theme(style="whitegrid")
    plt.figure(dpi=300) 
    plt.tight_layout()
    g = sns.catplot(
        data=df, kind="bar",
        x="Treatment Arm", y=yaxis, hue="Timepoint",
        palette="muted", height=6, legend=False
    )
    
    g.despine(left=True)
    g.set_ylabels(ylabel)
    #g.set_axis_labels("", "Body mass (g)")
    #plt.legend(loc='upper right')
    #g.legend.set_title("Avg. Mutations per Sample")
    
    
    #plt.figure(dpi=300)   
    plt.tight_layout()
    plt.legend(title="Timepoint", bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    g.savefig(ylabel + ".pdf")
    
    
###########################################################################3


def makeOncoprinterFiles(rowset, timepoint):
    
    # sort by sort by disease, subtype, treatment arm, progression, time till progression, 
    #        timepoint, mutation count
    
    df = rowset.copy()
    df["Time-till-Progression"] = pd.to_numeric(df["Time-till-Progression"])
    
    
    # fix variant classification
    # Oncoprinter: "MISSENSE", "INFRAME", "TRUNC", "SPLICE", "PROMOTER", or "OTHER" for a mutation alteration
    
    #df["Variant_Classification-Original"] = df["Variant_Classification"]
    df["Variant_Classification"] = "OTHER"
    
    '''
    df.loc[df["Variant_Classification-Original"] == "Missense_Mutation", 
                            "Variant_Classification"] = "MISSENSE"
    df.loc[df["Variant_Classification-Original"] == "In_Frame_Del", 
                            "Variant_Classification"] = "INFRAME"
    df.loc[df["Variant_Classification-Original"] == "In_Frame_Ins", 
                            "Variant_Classification-Original"] = "INFRAME"
    df.loc[df["Variant_Classification"] == "Frame_Shift_Del", 
                            "Variant_Classification-Original"] = "TRUNC"
    df.loc[df["Variant_Classification-Original"] == "Frame_Shift_Ins", 
                            "Variant_Classification"] = "TRUNC"
    df.loc[df["Variant_Classification-Original"] == "Nonsense_Mutation", 
                            "Variant_Classification"] = "TRUNC"
    df.loc[df["Variant_Classification-Original"] == "Splice_Site", 
                            "Variant_Classification"] = "SPLICE"
    df.loc[df["Variant_Classification-Original"] == "5'Flank", 
                            "Variant_Classification"] = "PROMOTER"
    df.loc[df["Variant_Classification-Original"] == "5'Flank",
                           "HGVSp_Short"] = "PROMOTER"
    '''
                            
    # take care of SVs that have verbose HGSVP short
    # if variant classification == TRA or HGSVP_short contains fusion, it's a fusion
    df.loc[(df["HGVSp_Short"].str.contains("fusion")) 
                     | (df["HGVSp_Short"].str.contains("Fusion")), 
                            ["Variant_Classification", "HGVSp_Short"]] = ["FUSION", "FUSION"]
    '''
    df.loc[df["Variant_Classification-Original"] == "TRA", 
                            ["Variant_Classification", "HGVSp_Short"]] = ["FUSION", "FUSION"]
    df.loc[df["Variant_Classification-Original"] == "DEL", 
                            ["Variant_Classification", "HGVSp_Short"]] = ["CNA", "HETLOSS"]
    '''                        
    
    df["HGVSp_Short"] = df["HGVSp_Short"].apply(lambda x: x.replace(" ", "_"))
    df["Timepoint"] = timepoint
    
    df['Unaltered'] = df.groupby('Sample')['AF'].transform('max')
    df['Unaltered'] = df['Unaltered'].apply(lambda x: "Yes" if x==0 else "No")
    unaltered = df.loc[df['Unaltered']=="Yes", ["Sample"]].drop_duplicates()
    m = df[(df["Unaltered"]=="No") & (df["AF"]>0)].copy()
    
    print(unaltered)
    
    # make and write mutation data
    # sort genes by number of occurrences
    m['geneFrequency'] = m.groupby('Hugo_Symbol')['Sample'].transform('nunique')
    m = m.sort_values('geneFrequency', ascending=False)
    mutationData = m[["Sample", "Hugo_Symbol", "HGVSp_Short", "Variant_Classification"]]
    mutationData = mutationData.append(unaltered)
    
    mutationData.to_csv("oncoprinter-mutations-" + timepoint + ".tsv", sep="\t", header=False, index=False)
    
    
    # make and write clinical data
    df = df.drop_duplicates(subset="Sample")
    #df["Subtype"] = df["Subtype"].apply(lambda x: x.replace(" ", "-"))
    df["Treatment Arm"] = df["Treatment Arm"].apply(lambda x: x.replace(" ", "-"))
    df["Timepoint"] = df["Timepoint"].apply(lambda x: x.replace(" ", "-"))
    
    
    # this doesn't really determine the oncoprinter sample order
    df = df.sort_values(["Disease", "Treatment Arm", "Progression",
                "Time-till-Progression", "Timepoint", "Sample"], 
                                             ascending = (True, True, True, False, True, True))
    df.rename(columns = {"Treatment Arm": "Treatment-Arm", "Time-till-Progression": "Time-till-Progression(number)"
                            }, inplace=True)
    
    df = df[["Sample", "Disease", "Treatment-Arm", "Progression", "Time-till-Progression(number)",
                "Timepoint",]]
    df.to_csv("oncoprinter-clinical-" + timepoint + ".tsv", sep="\t", index=False)
    
    
##########################################################################
    
def label_point(x, y, val, ax):
    a = pd.concat({'x': x, 'y': y, 'val': val}, axis=1)
    for i, point in a.iterrows():
        ax.text(point['x']+.02, point['y'], str(point['val']))
    
    

###########################################

def identifyEWPSamples(row, newrPatients):
    
    '''
    if the sample is a progression sample and the patient does not have an 8-week research sample and 
    the date of sample collection is between 6-12 weeks after randomizaton date, then the sample is 
    an "8-week/progression" sample 
    
    '''
    
    duration = (row['Research Lab Date'] - row['Randomization-Date']).days
    
    
    if (row["Timepoint"] == "Progression") and (row["Patient"] in newrPatients) and \
                (duration >= 42 and duration <= 84):
    
        return "8 weeks / Progression"
    
    else:
        
        return row["Timepoint"]
    



#############################################################################

def main():
    
    d = "/Users/patelj1/current/PromiseStudyJillian/analysis"
    mutationsFile = d + "/combined.tsv"
    
    os.chdir(d + "/breast")
    
    title = "Breast-SBRT"
    
    
    mutations = pd.read_csv(mutationsFile, sep="\t", index_col=False, keep_default_na=False, low_memory=False)
    
    '''
    # fix HGVSp_Short
    mutations.loc[(mutations["HGVSp_Short"]=="") | (mutations["HGVSp_Short"]=="-"), 
                  "HGVSp_Short"] = mutations["Variant_Classification"]
    '''
    
               
    # make Mutation field
    mutations["Mutation"] = mutations["Hugo_Symbol"] + "." + mutations["HGVSp_Short"]
    
    
    
    # fix AD, make 1 support 0, change AF accordingly
    mutations["AD"] = pd.to_numeric(mutations["AD"])
    mutations["DP"] = pd.to_numeric(mutations["DP"])
    mutations.loc[mutations["AD"]==1, "AD"] = 0
    mutations["AF"] = mutations["AD"]/mutations["DP"]
    
    # add Progression column
    mutations["Progression"] = mutations["Date of Progression"].apply(lambda x: "No" if x=="" else "Yes")
    
    '''
    # add needed columns
    #mutations["Timepoint-Original"] =  mutations["Timepoint"]
    #mutations.loc[mutations["Timepoint"]== "8 weeks / Progression", "Timepoint"] =  "Progression"
    
    mutations.loc[mutations["Disease"] == "Breast", "Subtype"] = "non-TNBC"
    mutations.loc[(mutations["Disease"]=="Breast") & ((mutations["ER"]=="-") &
                        (mutations["PR"]=="-") & (mutations["HER2"]=="-")), "Subtype"] = "TNBC"
                        
    mutations.loc[mutations["Disease"] == "Lung", "Subtype"] = "known driver"
    mutations.loc[(mutations["Disease"]=="Lung") & ((mutations["Marker"]=="Negative") |
                            (mutations["Marker"]=="SCC") | (mutations["Marker"]=="UNK")), "Subtype"] = "unknown driver"
    '''
    
    
    today = date.today().strftime("%m/%d/%y")
    mutations["Earliest Date of Progression"]  = mutations["Date of Progression"].str.strip()
    mutations.loc[mutations["Earliest Date of Progression"]=="", 
                                "Earliest Date of Progression"] = mutations["Off Study Date"].str.strip()
    mutations.loc[mutations["Earliest Date of Progression"]=="", 
                                "Earliest Date of Progression"] = today
                                
                                                
    mutations["Earliest Date of Progression"] = pd.to_datetime(mutations["Earliest Date of Progression"], errors='coerce')
    mutations["Randomization-Date"] = pd.to_datetime(mutations["Randomization Date"], errors='coerce')
    mutations["Research Lab Date"] = pd.to_datetime(mutations["Research Lab Date"], errors='coerce')
    
    mutations["Time-till-Progression"] = (mutations['Earliest Date of Progression'] - 
                                                mutations['Randomization-Date']).dt.days
    
    '''
     # mark 8-week/progression samples
    allPatients = mutations["Patient"].unique().tolist()
    ewrPatients = mutations.loc[mutations["Timepoint"]=="8W Research Labs", "Patient"].unique().tolist()
    newrPatients = np.setdiff1d(allPatients, ewrPatients)
    mutations["Timepoint"] = mutations.apply(
                                lambda row: identifyEWPSamples(row, newrPatients), axis=1)
    '''
    
    
    
    # filter mutations as needed
    mutations = mutations[(mutations["Disease"]=="Breast") & 
                         (mutations["Treatment Arm"]=="SBRT")]
    
    
    # save                                          
    mutations.to_csv("combined-df.tsv", sep="\t", index=False)
                        

   
    
    pd.options.display.max_columns = None
    pd.options.display.max_rows = None
    
    
    # get the different timepoints
    baseline = mutations.loc[(mutations["Timepoint"]=="Baseline")].copy()
    eightWeeks = mutations.loc[((mutations["Timepoint"]=="8 Weeks") |
                               (mutations["Timepoint"]=="8 Weeks-Progression"))].copy()
    progression = mutations.loc[((mutations["Timepoint"]=="Progression") |
                                 (mutations["Timepoint"]=="8 Weeks-Progression"))].copy()
    
    
    eightWeeks["Timepoint"] = "8 weeks"
    progression["Timepoint"] = "Progression"
    
    
    #makeOncoprinterFiles(baseline, "Baseline")
    #makeOncoprinterFiles(eightWeeks, "8-week")
    #makeOncoprinterFiles(progression, "Progression")
    
    
    # don't plot LOC101928728
    baseline = baseline[baseline["Hugo_Symbol"]!="LOC101928728"]
    eightWeeks = eightWeeks[eightWeeks["Hugo_Symbol"]!="LOC101928728"]
    progression = progression[progression["Hugo_Symbol"]!="LOC101928728"]
    
    
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
    
    # make df for trend plots
    # 2 types: baseline vs. 8-week and baseline vs. progression
    
    
    pairs = baseline.append(eightWeeks)
    pairs["NumTimepoints"] = pairs.groupby("Patient")["Timepoint"].transform('nunique')
    pairs = pairs[pairs["NumTimepoints"]==2]
    pairs = pairs.sort_values(by=["Time-till-Progression", "Patient", "Timepoint"], ascending=[False, True, False])
    pairs["P"] = "(" + pairs["Patient"] + ", " + pairs["Progression"] + ", " \
                    + pairs["Time-till-Progression"].astype(str) + ")"
    pairs = pairs.reset_index(drop=True)
    
    
    pairs.to_csv(title + "-bve.tsv", sep="\t", index=False)
    
    #sns.set_theme(style="whitegrid")
    plt.figure(dpi=300)
    plt.tight_layout()
    
    g = sns.FacetGrid(pairs, col="P", col_wrap=4, height=2, sharey=False)
    g.map_dataframe(sns.lineplot, x="Timepoint", y="AF",
                    hue="Mutation", palette="tab10", sort=False)
                    
    
    
    patients = pairs["Patient"].unique().tolist()
    print("BvE Patients:")
    print(patients)
    for i in range(0, len(g.axes)):
        ax = g.axes[i]
        patient = patients[i]
        t = pairs[(pairs["Patient"]==patient) & (pairs["Called"]=="Yes")]
        sns.scatterplot(data=t, x="Timepoint", y="AF", marker="s", color="0", ax=ax)
        ax.legend(fontsize=4, markerscale=0.05, bbox_to_anchor=(1.10, 1))
        
        
    g.fig.subplots_adjust(wspace=4)
    
    g.fig.subplots_adjust(top=0.85)
    g.fig.suptitle(title + "-bve")
    g.savefig(title + "-bve.pdf")
    
    
    
    
    
    
    pairs = baseline.append(progression)
    pairs["NumTimepoints"] = pairs.groupby("Patient")["Timepoint"].transform('nunique')
    pairs = pairs[pairs["NumTimepoints"]==2]
    pairs = pairs.sort_values(by=["Time-till-Progression", "Patient", "Timepoint"], ascending=[False, True, True])
    pairs["P"] = "(" + pairs["Patient"] + ", " + pairs["Progression"] + ", " + pairs["Time-till-Progression"].astype(str) + ")"
    pairs = pairs.reset_index(drop=True)
    
    
    pairs.to_csv(title + "-bvp.tsv", sep="\t", index=False)
    
    #sns.set_theme(style="whitegrid")
    plt.figure(dpi=300)
    plt.tight_layout()
    
    g = sns.FacetGrid(pairs, col="P", col_wrap=4, height=2, sharey=False)
    g.map_dataframe(sns.lineplot, x="Timepoint", y="AF",
                    hue="Mutation", palette="tab10", sort=False)
                    
    
    
    patients = pairs["Patient"].unique().tolist()
    print("BvP Patients:")
    print(patients)
    for i in range(0, len(g.axes)):
        ax = g.axes[i]
        patient = patients[i]
        t = pairs[(pairs["Patient"]==patient) & (pairs["Called"]=="Yes")]
        sns.scatterplot(data=t, x="Timepoint", y="AF", marker="s", color="0", ax=ax)
        ax.legend(fontsize=4, markerscale=0.05, bbox_to_anchor=(1.10, 1))
        
        
    g.fig.subplots_adjust(wspace=4)
    
    g.fig.subplots_adjust(top=0.85)
    g.fig.suptitle(title + "-bvp")
    g.savefig(title + "-bvp.pdf")
    
    
    
    
    
    
    









###################################################################################
   
   
    # get the sample lists
    #samples2 = patients.loc[(patients["Disease"] == "Breast") & (patients["CMO-ID"]!= ""), "CMO-ID"].tolist()
    
    
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
    
    




if __name__ == '__main__':
    main()
    
    
    
    
    
