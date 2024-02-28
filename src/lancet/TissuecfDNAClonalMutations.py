'''
Author: Juber Patel

'''

import os 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from datetime import date
from scipy import stats
import colorcet as cc
import sys


def main():

    d = "/Users/patelj1/current/PromiseStudyJillian/analysis/baseline-8-weeks-ctdna-fraction"

    os.chdir(d)
    
    tissueFile = "baseline-8-weeks-ctdna-fractions.tsv"
    cfDNAFile = "cleaned-cfdna.tsv"

    
    tissue = pd.read_csv(tissueFile, sep="\t", index_col=False, keep_default_na=False, low_memory=False)

    cfDNA = pd.read_csv(cfDNAFile, sep="\t", index_col=False, keep_default_na=False, low_memory=False)

    cfDNA["Mutation"] = cfDNA["Mutation"].apply(lambda x: x.replace(":p.", "."))


    patients = tissue["Patient"].unique()

    for patient in patients:

        m1 = set(tissue[tissue["Patient"]==patient]["Mutation"].unique())
        m2 = set(cfDNA[cfDNA["Patient"]==patient]["Mutation"].unique())

        if len(m2)==0:
            continue

        if m1.isdisjoint(m2):
            print(patient + " is disjoint")
            print(m1)
            print(m2)










if __name__ == '__main__':
    main()
