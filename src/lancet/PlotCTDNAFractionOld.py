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
from matplotlib.ticker import ScalarFormatter




def main():
    
    d = "/Users/patelj1/current/PromiseStudyJillian/analysis"
    mutationsFile = d + "/combined-ctDNAFraction.tsv"
    
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
    mutations["Progression"] = mutations["Date of Progression"].apply(lambda x: "No Progression" if x=="" else "Progression")
    
    
    
    today = date.today().strftime("%m/%d/%y")
    mutations["Earliest Date of Progression"]  = mutations["Date of Progression"].str.strip()
    mutations.loc[mutations["Earliest Date of Progression"]=="", 
                                "Earliest Date of Progression"] = mutations["Off Study Date"].str.strip()
    mutations.loc[mutations["Earliest Date of Progression"]=="", 
                                "Earliest Date of Progression"] = today
                                
                                                
    mutations["Earliest Date of Progression"] = pd.to_datetime(mutations["Earliest Date of Progression"], errors='coerce')
    mutations["Randomization-Date"] = pd.to_datetime(mutations["Randomization Date"], errors='coerce')
    
    mutations["Time-till-Progression"] = (mutations['Earliest Date of Progression'] - 
                                                mutations['Randomization-Date']).dt.days
    
    
    
    # save                                          
    mutations.to_csv("combined-df.tsv", sep="\t", index=False)
                        
   
    
    pd.options.display.max_columns = None
    pd.options.display.max_rows = None
    
    
    
    
    mutations = mutations.sort_values(by=["Disease", "Treatment Arm", "Progression", "Time-till-Progression",
                                            "Patient", "Timepoint"], 
                                      ascending=[True, False, True, False, True, True])
    
    
    
    
    # make a scatterplot of highest-baseline-af vs. ctDNA fraction
    s = mutations[mutations["computed_t_alt_cn_chosen"]!=""].copy()
    
    s.to_csv("s.tsv", sep="\t", index=False)
    
    '''
    s["Start_Position_tissue"] = s["Start_Position_tissue"].astype(float)
    s["Start_Position_tissue"] = s["Start_Position_tissue"].astype(int)
    s["Start_Position_tissue"] = s["Start_Position_tissue"].astype(int)
    s["M1"] = s["Chromosome_tissue"] + ":" + s["Start_Position_tissue"] + ":" + s["Reference_Allele_tissue"] + \
                               ":" + s["Tumor_Seq_Allele2_tissue"                         
    print(s["M1"])
    '''
    
    s["ctDNAFraction"] = s["ctDNAFraction"].astype(float)
    s.loc[s["Waltz_total_t_alt_count"]==1, "ctDNAFraction"] = 0.0
    s.loc[s["ctDNAFraction"]==0, "ctDNAFraction"] = 0.0005
    s.loc[s["AF"]==0, "AF"] = 0.0005
    print(len(s.index))
    
    
    sns.set_theme(style="whitegrid")
    
    plt.figure(dpi=300)
    plt.tight_layout()
    
    
    f, ax = plt.subplots(figsize=(5, 5))
    
    g = sns.scatterplot(data=s, x="AF", y="ctDNAFraction", color="red", alpha=0.5, ci=None, ax=ax)
    
    '''
    ax.set_xscale("log", base=2)
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.set_yscale("log", base=2)
    ax.yaxis.set_major_formatter(ScalarFormatter())
    '''
    
    plt.tight_layout()
    plt.savefig("ctDNAFraction-vs-HighestBaselineAF.pdf")
    
    


    
    


if __name__ == '__main__':
    main()
    
    
    
    
    