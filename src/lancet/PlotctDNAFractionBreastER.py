'''
Author: Juber Patel

'''

import os 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import date
from matplotlib import rcParams
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
from scipy import stats
import colorcet as cc
import sys


def wilcoxonTest(df, er, arm):

    
    t1 = df[(df["ER"]==er) & (df["TreatmentArm"]==arm) & (df["Timepoint"]=="Baseline")]["ctDNAFraction"]
    t2 = df[(df["ER"]==er) & (df["TreatmentArm"]==arm) & (df["Timepoint"]=="8 Weeks")]["ctDNAFraction"]

    print(er + ", " + arm + ": ")
    print(stats.wilcoxon(t1, t2))



def main():

    d = "/Users/patelj1/current/PromiseStudyJillian/analysis/baseline-8-weeks-ctdna-fraction"

    os.chdir(d)

    #ctDNAFractionFile = "ctDNA-fractions.tsv"

    cfDNAFractionsFile = "cfDNA-fractions.tsv"

    #ctDNAFractions = pd.read_csv(ctDNAFractionFile, sep="\t", index_col=False, keep_default_na=False, low_memory=False)

    erFile = "ER-PR-HER2.tsv"

    


    cfDNAFractions = pd.read_csv(cfDNAFractionsFile, sep="\t", index_col=False, keep_default_na=False, low_memory=False)

    er = pd.read_csv(erFile, sep="\t", index_col=False, keep_default_na=False, low_memory=False)


    cfDNAFractions = cfDNAFractions[cfDNAFractions["Disease"]=="Breast"]

    # rename the median AF column
    #cfDNAFractions.rename(columns={"MedianAF": "cfDNAMetric"}, inplace=True)

    # get the trend
    cfDNAFractions["T"] = cfDNAFractions["cfDNAMetric"]
    cfDNAFractions.loc[cfDNAFractions["Timepoint"] == "Baseline", "T"] = -cfDNAFractions["T"]
    cfDNAFractions["Trend"] = cfDNAFractions.groupby("Patient")["T"].transform('sum')
    cfDNAFractions.drop(columns=["T"], inplace=True)


    er = er[er["Disease"]=="Breast"]
    er = er[['Patient', 'ER', 'PR',	'HER2']]
    er.loc[er["ER"].str.contains("low"), "ER"] = "+"
    er = er.drop_duplicates(subset=['Patient'])

    cfDNAFractions = pd.merge(cfDNAFractions, er, on="Patient", how="left")

    cfDNAFractions.to_csv("x.tsv", sep="\t", index=False)
    

    # compute wilcoxon values
    cfDNAFractions["ctDNAFraction"] = cfDNAFractions["cfDNAMetric"]
    wilcoxonTest(cfDNAFractions, "+", "SBRT")
    wilcoxonTest(cfDNAFractions, "+", "No SBRT")
    wilcoxonTest(cfDNAFractions, "-", "SBRT")
    wilcoxonTest(cfDNAFractions, "-", "No SBRT")

    # drop the ctDNA fraction column
    cfDNAFractions.drop(columns=["ctDNAFraction"], inplace=True)
    

    #print(stats.ttest_rel(t1 ,t2))
    
    #print(stats.mannwhitneyu(t1, t2, method='asymptotic'))
    
    
    
    os.chdir("plots-ER")

    # make plots
    rcParams['figure.figsize'] = 5,6
    sns.set_theme(style="white")
    plt.figure(dpi=300)
    plt.tight_layout()

    

    
    # cfDNA-based measure
    for er in ["+", "-"]:
        for treatmentArm in ["SBRT", "No SBRT"]:
            df = cfDNAFractions.copy()
            df = df[(df["ER"]==er) & (df["TreatmentArm"]==treatmentArm)]
            #df = df[df["MedianAF"]<=0.02]
            numPatients = len(df["Patient"].unique())
            #p = sns.color_palette(cc.glasbey, n_colors=numPatients)
            title = "cfDNAMetric-ER" + er + "-" + treatmentArm
            df.loc[df["cfDNAMetric"]==0.0, "cfDNAMetric"] = 0.0001
            g = sns.lineplot(data = df, x="Timepoint", y="cfDNAMetric", hue="Patient")
            g.set(title = title)

            ax = plt.gca()
            ax.set_yscale("log", base=10)
            ax.yaxis.set_major_formatter(FormatStrFormatter('%.4f'))

            plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
            plt.tight_layout()
            plt.savefig(title + ".pdf")
            plt.close()

            # count the trends
            # remove cases where baseline and 8 weeks are both 0
            p1 = df[(df["Timepoint"]=="Baseline") & (df["cfDNAMetric"]>0.0)]["Patient"]
            p2 = df[(df["Timepoint"]=="8 Weeks") & (df["cfDNAMetric"]>0.0)]["Patient"]
            df = df[(df["Patient"].isin(p1)) | (df["Patient"].isin(p2))]
            up = len(df[df["Trend"]>0.0]["Patient"].unique())
            down = len(df[df["Trend"]<0.0]["Patient"].unique())
            print(er, treatmentArm, up, down, up+down)




    
    

    
            
            






if __name__ == "__main__":
    main()
