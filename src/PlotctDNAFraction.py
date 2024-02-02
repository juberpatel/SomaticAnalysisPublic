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


def wilcoxonTest(df, disease, arm):

    
    t1 = df[(df["Disease"]==disease) & (df["TreatmentArm"]==arm) & (df["Timepoint"]=="Baseline")]["ctDNAFraction"]
    t2 = df[(df["Disease"]==disease) & (df["TreatmentArm"]==arm) & (df["Timepoint"]=="8 Weeks")]["ctDNAFraction"]

    print(disease + ", " + arm + ": ")
    print(stats.wilcoxon(t1, t2))



def main():

    d = "/Users/patelj1/current/PromiseStudyJillian/analysis/baseline-8-weeks-ctdna-fraction"

    os.chdir(d)

    ctDNAFractionFile = "ctDNA-fractions.tsv"

    cfDNAFractionsFile = "cfDNA-fractions.tsv"

    ctDNAFractions = pd.read_csv(ctDNAFractionFile, sep="\t", index_col=False, keep_default_na=False, low_memory=False)


    cfDNAFractions = pd.read_csv(cfDNAFractionsFile, sep="\t", index_col=False, keep_default_na=False, low_memory=False)

    # rename the median AF column
    cfDNAFractions.rename(columns={"MedianAF": "cfDNAMetric"}, inplace=True)

    
    # get the trend
    cfDNAFractions["T"] = cfDNAFractions["cfDNAMetric"]
    cfDNAFractions.loc[cfDNAFractions["Timepoint"] == "Baseline", "T"] = -cfDNAFractions["T"]
    cfDNAFractions["Trend"] = cfDNAFractions.groupby("Patient")["T"].transform('sum')
    cfDNAFractions.drop(columns=["T"], inplace=True)
     
    

    
    # disease, treatment arm, 
    # compute wilcoxon values

    wilcoxonTest(ctDNAFractions, "Lung", "SBRT")
    wilcoxonTest(ctDNAFractions, "Lung", "No SBRT")
    wilcoxonTest(ctDNAFractions, "Breast", "SBRT")
    wilcoxonTest(ctDNAFractions, "Breast", "No SBRT")

    cfDNAFractions["ctDNAFraction"] = cfDNAFractions["cfDNAMetric"]
    wilcoxonTest(cfDNAFractions, "Lung", "SBRT")
    wilcoxonTest(cfDNAFractions, "Lung", "No SBRT")
    wilcoxonTest(cfDNAFractions, "Breast", "SBRT")
    wilcoxonTest(cfDNAFractions, "Breast", "No SBRT")

    # drop the ctDNA fraction column
    cfDNAFractions.drop(columns=["ctDNAFraction"], inplace=True)



    '''
    t1 = ctDNAFractions[(ctDNAFractions["Disease"]=="Lung") & (ctDNAFractions["TreatmentArm"]=="SBRT") & (ctDNAFractions["Timepoint"]=="Baseline")]["ctDNAFraction"]

    t2 = ctDNAFractions[(ctDNAFractions["Disease"]=="Lung") & (ctDNAFractions["TreatmentArm"]=="SBRT") & (ctDNAFractions["Timepoint"]=="8 Weeks")]["ctDNAFraction"]

    
    t1 = cfDNAFractions[(cfDNAFractions["Disease"]=="Lung") & (cfDNAFractions["TreatmentArm"]=="SBRT") & (cfDNAFractions["Timepoint"]=="Baseline")]["MedianAF"]

    t2 = cfDNAFractions[(cfDNAFractions["Disease"]=="Lung") & (cfDNAFractions["TreatmentArm"]=="SBRT") & (cfDNAFractions["Timepoint"]=="8 Weeks")]["MedianAF"]
    '''

    

    #print(stats.ttest_rel(t1 ,t2))
    
    #print(stats.mannwhitneyu(t1, t2, method='asymptotic'))
    
    
    
    os.chdir("plots")

    # make plots
    rcParams['figure.figsize'] = 5,6
    sns.set_theme(style="white")
    plt.figure(dpi=300)
    plt.tight_layout()

    '''
    # ctDNA fraction
    for disease in ["Lung", "Breast"]:
        for treatmentArm in ["SBRT", "No SBRT"]:
            df = ctDNAFractions[(ctDNAFractions["Disease"]==disease) & (ctDNAFractions["TreatmentArm"]==treatmentArm)]
            #df = df[df["ctDNAFraction"]<=0.04]
            numPatients = len(df["Patient"].unique())
            #p = sns.color_palette(cc.glasbey, n_colors=numPatients)
            title = "ctDNAFraction-" + disease + "-" + treatmentArm
            df.loc[df["ctDNAFraction"]==0.0, "ctDNAFraction"] = 0.0001
            g = sns.lineplot(data = df, x="Timepoint", y="ctDNAFraction", hue="Patient")
            g.set(title = title)

            ax = plt.gca()
            ax.set_yscale("log", base=10)
            ax.yaxis.set_major_formatter(FormatStrFormatter('%.4f'))

            plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
            plt.tight_layout()
            plt.savefig(title + ".pdf")
            plt.close()

    
    # cfDNA-based measure
    for disease in ["Lung", "Breast"]:
        for treatmentArm in ["SBRT", "No SBRT"]:
            df = cfDNAFractions.copy()
            df = df[(df["Disease"]==disease) & (df["TreatmentArm"]==treatmentArm)]
            #df = df[df["MedianAF"]<=0.02]
            numPatients = len(df["Patient"].unique())
            #p = sns.color_palette(cc.glasbey, n_colors=numPatients)
            title = "cfDNAMetric-" + disease + "-" + treatmentArm
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
            print(disease, treatmentArm, up, down, up+down)

    '''


    


    # comparison of ctDNA fraction and median AF

    
    t = cfDNAFractions[['Patient', 'Timepoint', 'cfDNAMetric']]
    ctDNAFractions = ctDNAFractions.merge(t, on=['Patient', 'Timepoint'], how='left')
    ctDNAFractions["fc"] = 0.0
    ctDNAFractions.loc[ctDNAFractions["cfDNAMetric"]!=0.0, "fc"] = ctDNAFractions["ctDNAFraction"]/ctDNAFractions["cfDNAMetric"]
    ctDNAFractions["diff"] = ctDNAFractions["ctDNAFraction"] - ctDNAFractions["cfDNAMetric"]
    ctDNAFractions = ctDNAFractions[(ctDNAFractions["Timepoint"]=="Baseline") | (~ctDNAFractions["Patient"].isin(["P083", "P041", "P067"]))]
    ctDNAFractions.to_csv("t.tsv", sep="\t", index=False)

    
    # scatterplot
    rcParams['figure.figsize'] = 5,5
    sns.set_theme(style="whitegrid")
    title = "ctDNAFraction-cfDNAMetric-scatterplot"

    df = ctDNAFractions.copy()
    
    
    #df.loc[df["ctDNAFraction"]==0.0, "ctDNAFraction"] = 0.0001
    #df.loc[df["cfDNAMetric"]==0.0, "cfDNAMetric"] = 0.0001
    g = sns.scatterplot(data = df, x="ctDNAFraction", y="cfDNAMetric")
    g.set(title = title)

    ax = plt.gca()
    ax.set_xlim(0, 0.4)
    ax.set_ylim(0, 0.8)
    #ax.set_xscale("log", base=10)
    #ax.set_yscale("log", base=10)
    #ax.xaxis.set_major_formatter(FormatStrFormatter('%.4f')) #(ScalarFormatter('.4f'))
    #ax.yaxis.set_major_formatter(FormatStrFormatter('%.4f'))


    #plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
    plt.tight_layout()
    plt.savefig(title + ".pdf")
    plt.close()
    

    # Pearson correlation
    print(stats.pearsonr(df["ctDNAFraction"], df["cfDNAMetric"]))
    

    ctDNAFractions = ctDNAFractions.melt(id_vars=['Patient', 'Timepoint', 'Disease', 'TreatmentArm',], value_vars=['ctDNAFraction', 'cfDNAMetric'], var_name='Method', value_name='Value')

    '''
    # timepoint plots
    rcParams['figure.figsize'] = 5,6
    sns.set_theme(style="white")
    for disease in ["Lung", "Breast"]:
        for treatmentArm in ["SBRT", "No SBRT"]:
            df = ctDNAFractions.copy()
            df = df[(df["Disease"]==disease) & (df["TreatmentArm"]==treatmentArm)]
            #df = df[df["Value"]<=0.05]
            df.loc[df["Value"]==0.0, "Value"] = 0.0001
            numPatients = len(df["Patient"].unique())
            #p = sns.color_palette(cc.glasbey, n_colors=numPatients)

            title = "ctDNAFraction-cfDNAMetric" + disease + "-" + treatmentArm
            g = sns.lineplot(data = df, x="Timepoint", y="Value", hue="Patient", style="Method")
            g.set(title = title)

            
            ax = plt.gca()
            ax.set_yscale("log", base=10)
            ax.yaxis.set_major_formatter(FormatStrFormatter('%.4f')) #(ScalarFormatter('.4f'))
            
            plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
            plt.tight_layout()
            plt.savefig(title + ".pdf")
            plt.close()
    '''
            
    

    






if __name__ == "__main__":
    main()