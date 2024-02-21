'''
Created on Jul 15, 2021

@author: patelj1
'''

import os 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import date
from matplotlib import rcParams
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter




def makeOncoprinterFiles(rowset, timepoint, disease):
    
    # sort by sort by disease, subtype, treatment arm, progression, time till progression, 
    #        timepoint, mutation count
    
    df = rowset.copy()
    df["PFS"] = pd.to_numeric(df["PFS"])
    
    
    # fix variant classification
    # Oncoprinter: "MISSENSE", "INFRAME", "TRUNC", "SPLICE", "PROMOTER", or "OTHER" for a mutation alteration
    
    # make silent mutations' AF 0
    df.loc[df["Variant_Classification"]=="Silent", "AF"] = 0
    
    df["Variant_Classification-Original"] = df["Variant_Classification"]
    df["Variant_Classification"] = "OTHER"
    
    
    df.loc[df["Variant_Classification-Original"] == "Missense_Mutation", 
                            "Variant_Classification"] = "MISSENSE"
    df.loc[df["Variant_Classification-Original"] == "In_Frame_Del", 
                            "Variant_Classification"] = "INFRAME"
    df.loc[df["Variant_Classification-Original"] == "In_Frame_Ins", 
                            "Variant_Classification"] = "INFRAME"
    df.loc[df["Variant_Classification-Original"] == "Frame_Shift_Del", 
                            "Variant_Classification"] = "TRUNC"
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
    
                            
    # take care of SVs that have verbose HGSVP short
    # if variant classification == TRA or HGSVP_short contains fusion, it's a fusion
    df.loc[(df["HGVSp_Short"].str.contains("fusion")) 
                     | (df["HGVSp_Short"].str.contains("Fusion")), 
                            ["Variant_Classification", "HGVSp_Short"]] = ["FUSION", "FUSION"]
    
    df.loc[df["Variant_Classification-Original"] == "TRA", 
                            ["Variant_Classification", "HGVSp_Short"]] = ["FUSION", "FUSION"]
    df.loc[df["Variant_Classification-Original"] == "DEL", 
                            ["Variant_Classification", "HGVSp_Short"]] = ["CNA", "HETLOSS"]
    
    
    df["HGVSp_Short"] = df["HGVSp_Short"].apply(lambda x: x.replace(" ", "_"))
    df["Timepoint"] = timepoint
    df["Disease"] = disease
    
    df['Unaltered'] = df.groupby('Sample')['AF'].transform('max')
    df['Unaltered'] = df['Unaltered'].apply(lambda x: "Yes" if x==0 else "No")
    unaltered = df.loc[df['Unaltered']=="Yes", ["Sample"]].drop_duplicates()
    m = df[(df["Unaltered"]=="No") & (df["AF"]>0)].copy()
    
    #print(unaltered)
    
    # make and write mutation data
    # sort genes by number of occurrences
    m['geneFrequency'] = m.groupby('Hugo_Symbol')['Sample'].transform('nunique')
    #m = m[m['geneFrequency']>=2]
    m = m.sort_values('geneFrequency', ascending=False)
    mutationData = m[["Sample", "Hugo_Symbol", "HGVSp_Short", "Variant_Classification"]]
    mutationData = mutationData.append(unaltered)
    
    mutationData.to_csv("oncoprinter-mutations-" + timepoint + "-" + disease + ".tsv", sep="\t", header=False, index=False)
    
    
    # make and write clinical data
    
    #calculate number of mutations per sample or per patient-timepoint, after removing silent mutations
    df = df[df["Variant_Classification"]!="Silent"]
    df["MutationCount"] = df.groupby("Sample").M.transform('count')
    
    # remove 
    
    
    df = df.drop_duplicates(subset="Sample")
    #df["Subtype"] = df["Subtype"].apply(lambda x: x.replace(" ", "-"))
    df["Treatment Arm"] = df["Treatment Arm"].apply(lambda x: x.replace(" ", "-"))
    df["Progression"] = df["Progression"].apply(lambda x: x.replace(" ", "-"))
    #df["Timepoint"] = df["Timepoint"].apply(lambda x: x.replace(" ", "-"))
    
    
    # this doesn't really determine the oncoprinter sample order
    df = df.sort_values(["Treatment Arm", "Progression",
                "PFS", "MutationCount", "Sample"], 
                                             ascending = (True, True, False, True, True))
    df.rename(columns = {"Treatment Arm": "Treatment-Arm", "PFS": "PFS(number)", "MutationCount": "MutationCount(number)"
                            }, inplace=True)
    
    df = df[["Sample", "Treatment-Arm", "Progression", "PFS(number)",
                "MutationCount(number)"]]
    df.to_csv("oncoprinter-clinical-" + timepoint + "-" + disease + ".tsv", sep="\t", index=False)
    
    
##########################################################################
    





def main():
    
    d = "/Users/patelj1/current/PromiseStudyJillian/analysis"
    mutationsFile = d + "/merged-AFGM.tsv"
    
    os.chdir(d + "/summary")
    
    title = "Summary"
    
    
    mutations = pd.read_csv(mutationsFile, sep="\t", index_col=False, keep_default_na=False, low_memory=False)
    
    # make Mutation field
    mutations["Mutation"] = mutations["Hugo_Symbol"] + "." + mutations["HGVSp_Short"]
    
    
    
    # fix AD, make 1 support 0, change AF accordingly
    mutations["AD"] = pd.to_numeric(mutations["AD"])
    mutations["DP"] = pd.to_numeric(mutations["DP"])
    mutations.loc[mutations["AD"]==1, "AD"] = 0
    mutations["AF"] = mutations["AD"]/mutations["DP"]
    
    # add Progression column
    mutations["Progression"] = mutations["POD"].apply(lambda x: "Progression" if x=="Y" else "No Progression")
    
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
    
    '''
    # change this accordingly
    #today = date.today().strftime("%m/%d/%y")
    today = "08/02/21"
    mutations["Earliest Date of Progression"]  = mutations["Date of Progression"].str.strip()
    mutations.loc[mutations["Earliest Date of Progression"]=="", 
                                "Earliest Date of Progression"] = mutations["Off Study Date"].str.strip()
    mutations.loc[mutations["Earliest Date of Progression"]=="", 
                                "Earliest Date of Progression"] = today
                                
                                                
    mutations["Earliest Date of Progression"] = pd.to_datetime(mutations["Earliest Date of Progression"], errors='coerce')
    mutations["Randomization-Date"] = pd.to_datetime(mutations["Randomization Date"], errors='coerce')
    
    mutations["Time-till-Progression"] = (mutations['Earliest Date of Progression'] - 
                                                mutations['Randomization-Date']).dt.days
    '''
    
    
    mutations['PFSDate'] = pd.to_datetime(mutations["PFSDate"], errors='coerce')
    mutations['RandomizationDate'] = pd.to_datetime(mutations["RandomizationDate"], errors='coerce')
    
    mutations["PFS"] = (mutations['PFSDate'] - mutations['RandomizationDate']).dt.days
    
    
    # save                                          
    #mutations.to_csv("combined-df.tsv", sep="\t", index=False)
                        
   
    
    pd.options.display.max_columns = None
    pd.options.display.max_rows = None
    
    # get the different timepoints
    baseline = mutations.loc[(mutations["Timepoint"]=="Baseline")].copy()
    #eightWeeks = mutations.loc[((mutations["Timepoint"]=="8 weeks") |
    #                           (mutations["Timepoint"]=="8 weeks / Progression"))].copy()
    progression = mutations.loc[((mutations["Timepoint"]=="Progression") |
                                 (mutations["Timepoint"]=="8 weeks / Progression"))].copy()
    
    
    #eightWeeks["Timepoint"] = "8 weeks"
    progression["Timepoint"] = "Progression"
    
    
    lung = baseline[baseline["Disease"]=="Lung"]
    breast = baseline[baseline["Disease"]=="Breast"]
    makeOncoprinterFiles(lung, "Baseline", "Lung")
    makeOncoprinterFiles(breast, "Baseline", "Breast")
    
    #makeOncoprinterFiles(eightWeeks, "8-week")
    #makeOncoprinterFiles(progression, "Progression")
    
    
    # don't plot LOC101928728
    # mutations = mutations[mutations["Hugo_Symbol"]!="LOC101928728"]
    # baseline = baseline[baseline["Hugo_Symbol"]!="LOC101928728"]
    #eightWeeks = eightWeeks[eightWeeks["Hugo_Symbol"]!="LOC101928728"]
    # progression = progression[progression["Hugo_Symbol"]!="LOC101928728"]
    
    
    
    
    
    mutations["Dummy"] = 1.0
    #mutations.loc[mutations["Progression"]=="Yes", "Progression"] = "Progression"
    #mutations.loc[mutations["Progression"]=="No", "Progression"] = "No Progression"
    
    mutations = mutations.sort_values(by=["Disease", "Treatment Arm", "Progression", "PFS",
                                            "Patient", "Timepoint"], 
                                      ascending=[True, False, True, False, True, True])
    
    
    mutations.to_csv("Mean-AF.tsv", sep="\t", index=False)
    
    
    
    sns.set_theme(style="whitegrid")
    plt.figure(dpi=300)
    plt.tight_layout()
    
    fig, axs = plt.subplots(6,1,
                      figsize=(50,25),
                      sharex=True,
                      gridspec_kw=dict(height_ratios=[1, 1, 1, 4, 4, 4]))
    
    
    #0
    palette = {"Breast":"crimson",
               "Lung":"green"}
    
    g = sns.barplot(
        data=mutations,
        x="Patient", y="Dummy", hue="Disease", dodge=False, ci=None,
        palette=palette, ax=axs[0]
    )
    
    g.set(yticklabels = [])
    g.set(ylabel=None)
    g.set(xlabel=None)
    
    #1
    palette = {"SBRT":"crimson",
               "No SBRT":"green"}
    
    g = sns.barplot(
        data=mutations,
        x="Patient", y="Dummy", hue="Treatment Arm", dodge=False, ci=None,
        palette=palette, ax=axs[1]
    )
    
    g.set(yticklabels = [])
    g.set(ylabel=None)
    g.set(xlabel=None)
    
    
    #2
    palette = {"Progression":"crimson",
               "No Progression":"green"}
    
    g = sns.barplot(
        data=mutations,
        x="Patient", y="Dummy", hue="Progression", dodge=False, ci=None,
        palette=palette, ax=axs[2]
    )
    
    g.set(yticklabels = [])
    g.set(ylabel=None)
    g.set(xlabel=None)
    
    
    #3
    g = sns.barplot(
        data=mutations,
        x="Patient", y="PFS", dodge=False, ci=None,
        color="darkorchid", ax=axs[3]
    )
    
    g.set_ylabel(ylabel="PFS", size=30)
    g.tick_params(labelsize=30)
    g.set(xlabel=None)
    
    
    
    #4
    g = sns.barplot(
        data=mutations,
        x="Patient", y="AFMedian", hue="Timepoint",
        palette="bright", ax=axs[4] #aspect=40/5, 
    )
    
    axs[4].set_yscale("log", base=2)
    axs[4].yaxis.set_major_formatter(ScalarFormatter())
    
    g.set_ylabel(ylabel="Median AF", size=30)
    g.tick_params(labelsize=30)
    g.set(xlabel=None)
    
    
    
    #5
    g = sns.barplot(
        data=mutations,
        x="Patient", y="AFGM", hue="Timepoint",
        palette="bright", ax=axs[5] 
    )
    
    axs[5].set_yscale("log", base=2)
    axs[5].yaxis.set_major_formatter(ScalarFormatter())
    
    g.set_ylabel(ylabel="AF Geometric Mean", size=30)
    g.tick_params(labelsize=30)
    g.set_xlabel(xlabel="Patient", size=30)
    
    plt.xticks(size=20, rotation=90)
    #plt.yticks(size=30)

    
    # add legends
    axs[0].legend(fontsize=20, bbox_to_anchor=(-0.02, 1))
    axs[1].legend(fontsize=20, bbox_to_anchor=(-0.02, 1))
    axs[2].legend(fontsize=20, bbox_to_anchor=(-0.02, 1))
    #axs[3].legend(fontsize=20, bbox_to_anchor=(-0.02, 1))
    axs[4].legend(fontsize=20, bbox_to_anchor=(-0.02, 1))
    axs[5].legend(fontsize=20, bbox_to_anchor=(-0.02, 1))
    
    
    plt.tight_layout()
    plt.savefig("Mean-AF.pdf")
    

    
    
    # make boxplots of medain AF
    
    for disease in ["Lung", "Breast"]:
        for metric in ["AFGM", "AFMedian"]:
            
            # bve
            df = mutations.copy()
            df = df[df["Patient"]!="P065"]
            df = df[df["Disease"]==disease]
            df.loc[df["Timepoint"]=="8 Weeks-Progression", "Timepoint"] = "8 Weeks"
            df = df[df["Timepoint"]!="Progression"]
            df = df.drop_duplicates(subset=['Sample'])
            
            df = df[df.groupby('Patient').Sample.transform('count')>1].copy()
            
            
            df.loc[df[metric]==0, metric] = 0.0002
            
            df.to_csv(disease + "-" + metric + "-bve.tsv", sep="\t", index=False)
            
            sns.set_theme(style="whitegrid")
            plt.figure(dpi=300)
            plt.tight_layout()
            
            sns.boxplot(x="Treatment Arm", y=metric, hue="Timepoint", hue_order=["Baseline", "8 Weeks"], data=df,
                        color="white") #, palette="rainbow")
            g = sns.swarmplot(x='Treatment Arm', y=metric, hue='Timepoint', hue_order=["Baseline", "8 Weeks"], data=df, 
                              dodge=True, palette='viridis')
            g.set(title = disease + " " + metric + " across Treatment Arms")
            ax = plt.gca()
            ax.set_yscale("log", base=2)
            ax.yaxis.set_major_formatter(FormatStrFormatter('%.4f')) #(ScalarFormatter('.4f'))
            plt.tight_layout()
            plt.savefig(disease + "-" + metric + "-bve.pdf")
            plt.close()
            
            
            # bve-fc
            df = df.pivot_table(index=['Patient', 'Treatment Arm', 'PFS'], columns=['Timepoint'], 
                                values=metric)
            df.reset_index(inplace=True)
            df.columns.name = None
            df['fc'] = df['8 Weeks']/df['Baseline']
            print(df)
            df.to_csv(disease + "-" + metric + "-bve-fc.tsv", sep="\t", index=False)
            
            sns.set_theme(style="whitegrid")
            plt.figure(dpi=300)
            plt.tight_layout()
            
            sns.boxplot(x="Treatment Arm", y='fc', data=df,
                        color="white") #, palette="rainbow")
            g = sns.swarmplot(x='Treatment Arm', y='fc', data=df, 
                              dodge=True, palette='viridis')
            
            g.set(title=disease + " " + metric + " fold change across Treatment Arms")
            ax = plt.gca()
            ax.set_yscale("log", base=2)
            ax.yaxis.set_major_formatter(FormatStrFormatter('%.4f')) #(ScalarFormatter('.4f'))
            plt.tight_layout()
            plt.savefig(disease + "-" + metric + "-bve-fc.pdf")
            plt.close()
            
            
            # bve fc-PFS scatterplot
            sns.set_theme(style="whitegrid")
            plt.figure(dpi=300)
            plt.tight_layout()
            
            g = sns.scatterplot(data=df, x="fc", y="PFS", hue="Treatment Arm")
            
            g.set(title=disease + " " + metric + " fold change vs PFS")
            ax = plt.gca()
            ax.set_xscale("log", base=2)
            ax.xaxis.set_major_formatter(FormatStrFormatter('%.4f')) #(ScalarFormatter('.4f'))
            plt.tight_layout()
            plt.savefig(disease + "-" + metric + "-bve-fc-vs-pfs.pdf")
            plt.close()
            
            
            
            
            # bvp
            df = mutations.copy()
            df = df[df["Patient"]!="P065"]
            df = df[df["Disease"]==disease]
            df.loc[df["Timepoint"]=="8 Weeks-Progression", "Timepoint"] = "Progression"
            df = df[df["Timepoint"]!="8 Weeks"]
            df = df.drop_duplicates(subset=['Sample'])
            
            df = df[df.groupby('Patient').Sample.transform('count')>1].copy() 
            
            df.loc[df[metric]==0, metric] = 0.0002
            
            df.to_csv(disease + "-" + metric + "-bvp.tsv", sep="\t", index=False)
            
            sns.set_theme(style="whitegrid")
            plt.figure(dpi=300)
            plt.tight_layout()
            
            sns.boxplot(x="Treatment Arm", y=metric, hue="Timepoint", hue_order=["Baseline", "Progression"], data=df,
                        color="white") #, palette="rainbow")
            g = sns.swarmplot(x='Treatment Arm', y=metric, hue='Timepoint', hue_order=["Baseline", "Progression"], data=df, 
                              dodge=True, palette='viridis')
            g.set(title= disease + " " + metric + " across Treatment Arms")
            ax = plt.gca()
            ax.set_yscale("log", base=2)
            ax.yaxis.set_major_formatter(FormatStrFormatter('%.4f')) #(ScalarFormatter('.4f'))
            plt.tight_layout()
            plt.savefig(disease + "-" + metric + "-bvp.pdf")
            plt.close()
            
            
            # bvp-fc
            df = df.pivot_table(index=['Patient', 'Treatment Arm', 'PFS'], columns=['Timepoint'], 
                                values=metric)
            df.reset_index(inplace=True)
            df.columns.name = None
            df['fc'] = df['Progression']/df['Baseline']
            print(df)
            df.to_csv(disease + "-" + metric + "-bvp-fc.tsv", sep="\t", index=False)
            
            sns.set_theme(style="whitegrid")
            plt.figure(dpi=300)
            plt.tight_layout()
            
            sns.boxplot(x="Treatment Arm", y='fc', data=df,
                        color="white") #, palette="rainbow")
            g = sns.swarmplot(x='Treatment Arm', y='fc', data=df, 
                              dodge=True, palette='viridis')
            g.set(title=disease + " " + metric + " fold change across Treatment Arms")
            ax = plt.gca()
            ax.set_yscale("log", base=2)
            ax.yaxis.set_major_formatter(FormatStrFormatter('%.4f')) #(ScalarFormatter('.4f'))
            plt.tight_layout()
            plt.savefig(disease + "-" + metric + "-bvp-fc.pdf")
            plt.close()
            
            
            # bvp fc-PFS scatterplot
            sns.set_theme(style="whitegrid")
            plt.figure(dpi=300)
            plt.tight_layout()
            
            g = sns.scatterplot(data=df, x="fc", y="PFS", hue="Treatment Arm")
            
            g.set(title=disease + " " + metric + " fold change vs PFS")
            ax = plt.gca()
            ax.set_xscale("log", base=2)
            ax.xaxis.set_major_formatter(FormatStrFormatter('%.4f')) #(ScalarFormatter('.4f'))
            plt.tight_layout()
            plt.savefig(disease + "-" + metric + "-bvp-fc-vs-pfs.pdf")
            plt.close()
            
            
            
            
    


if __name__ == '__main__':
    main()
    
    
    
    
    