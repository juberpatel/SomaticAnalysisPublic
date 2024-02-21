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




def main():
    
    d = "/Users/patelj1/current/PromiseStudyJillian/analysis/baseline-8-weeks-ctdna-fraction"

    os.chdir(d)

    #ctDNAFractionFile = "ctDNA-fractions.tsv"

    cfDNAFractionsFile = "cfDNA-fractions.tsv"

    #ctDNAFractions = pd.read_csv(ctDNAFractionFile, sep="\t", index_col=False, keep_default_na=False, low_memory=False)
    

    cfDNAFractions = pd.read_csv(cfDNAFractionsFile, sep="\t", index_col=False, keep_default_na=False, low_memory=False)


    df = cfDNAFractions[(cfDNAFractions["Disease"]=="Lung") & (cfDNAFractions["TreatmentArm"]=="No SBRT") & (cfDNAFractions["Timepoint"]=="Baseline")]

    min = df["cfDNAMetric"].min()
    max = df["cfDNAMetric"].max()
    mean = df["cfDNAMetric"].mean()
    median = df["cfDNAMetric"].median()
    print("Lung No SBRT Baseline")
    print(min, max, mean, median)


    df = cfDNAFractions[(cfDNAFractions["Disease"]=="Lung") & (cfDNAFractions["TreatmentArm"]=="No SBRT") & (cfDNAFractions["Timepoint"]=="8 Weeks")]

    min = df["cfDNAMetric"].min()
    max = df["cfDNAMetric"].max()
    mean = df["cfDNAMetric"].mean()
    median = df["cfDNAMetric"].median()
    print("Lung No SBRT 8 Weeks")
    print(min, max, mean, median)

    
























if __name__ == "__main__":
    main()