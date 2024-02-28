'''
Created on Jul 28, 2021

@author: Juber Patel
'''

import os
import pandas as pd
from scipy import stats
import statistics



def computeGeometricMeanOld(s):
    
    nonZero = s[s!=0]
    
    if nonZero.count() == 0:
        return 0
    else:
        return stats.gmean(nonZero)





def computeGeometricMean(s):
    
    #print(s)
    
    additive = 1.0
    
    s1 = s + additive
    
    #print(s1)
    
    return stats.gmean(s1) - additive





def computeMedian(s):
    return statistics.median(s)


######################################################################


def main():
    
    d = "/Users/patelj1/current/PromiseStudyJillian/analysis"
    
    os.chdir(d)
    
    
    cfDNAFile = "merged.tsv"
    
    
    cfDNA = pd.read_csv(cfDNAFile, sep="\t", index_col=False, keep_default_na=False, low_memory=False)
    #cfDNA["M"] = cfDNA["Chromosome"].astype(str) + ":" + cfDNA["Start_Position"].astype(str) + ":" \
    #                         + cfDNA["Reference_Allele"] + ":" + cfDNA["Tumor_Seq_Allele2"]
    # drop old TUMOR_MAF, we need this column
    #cfDNA = cfDNA.drop("TUMOR_MAF", axis=1)
    
    
     
    # find mutations that have 0 AF across all samples in the patient
    removal = []
    g = cfDNA.groupby(["Patient", "M"])
    for name, group in g:
        m = group.AF.max()
        if m < 0.00001:
            removal.append(name)
            print(name, m)
            
    
    # remove mutations that have 0 AF across all samples in the patient
    filtered = cfDNA[~(cfDNA[['Patient', 'M']].apply(tuple, axis=1).isin(removal))]
    
    #print(filtered)
    
    
    # get geometric mean of AF
    filtered["AFGM"] = filtered.groupby(["Sample"]).AF.transform(computeGeometricMean)
                    
    # get median of AF
    filtered["AFMedian"] = filtered.groupby(["Sample"]).AF.transform(computeMedian)
    
    
    
    # save
    filtered.to_csv("merged-AFGM.tsv", sep="\t", index=False)
    
    
    
    
    
    
    
    
    
    






















if __name__ == '__main__':
    main()
