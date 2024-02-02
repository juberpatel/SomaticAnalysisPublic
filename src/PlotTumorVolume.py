import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import date
from matplotlib import rcParams
from matplotlib.ticker import ScalarFormatter
from pandas.api.types import CategoricalDtype
import sys
import colorcet as cc



def main():
    
    d = "/Users/patelj1/current/Beryl-17-180/analysis-new/tumor-volume-correlation/"

    os.chdir(d)

    ctDNAFile = "../seaborn-plots/all-with-ctDNAFraction-oncokb-processed.tsv"
    tumorVolumeFile = "tumor-volume-data.tsv"
    c1d1File = "../c1d1.tsv"

    ctDNA = pd.read_csv(ctDNAFile, sep="\t", index_col=False, keep_default_na=False, low_memory=False)
    tumorVolume = pd.read_csv(tumorVolumeFile, sep="\t", index_col=False, keep_default_na=False, low_memory=False)
    c1d1 = pd.read_csv(c1d1File, sep="\t", index_col=False, keep_default_na=False, low_memory=False)


    ctDNA = ctDNA[["Patient", "Patient1", "Timepoint", "ctDNAFraction"]]
    #ctDNA["ctDNAFraction"] = ctDNA["ctDNAFraction"].astype(str)
    ctDNA.drop_duplicates(subset=["Patient", "Timepoint"], inplace=True)


    c1d1 = c1d1[["Patient", "C1D1"]]
    c1d1.drop_duplicates(subset=["Patient"], inplace=True)
    c1d1["C1D1"] = pd.to_datetime(c1d1["C1D1"])


    tumorVolume = tumorVolume.merge(ctDNA, how="left", on=["Patient", "Timepoint"])
    tumorVolume = tumorVolume.merge(c1d1, how="left", on=["Patient"])
    tumorVolume.to_csv("t.tsv", sep="\t", index=False)

    
    tumorVolume["Volume (cm3)"] = pd.to_numeric(tumorVolume["Volume (cm3)"], errors="coerce")
    tumorVolume["ctDNAFraction"] = pd.to_numeric(tumorVolume["ctDNAFraction"], errors="coerce")
    tumorVolume = tumorVolume[(~tumorVolume["ctDNAFraction"].isnull()) | (~tumorVolume["Volume (cm3)"].isnull())]


    # compute days since first blood draw

    tumorVolume["cfDNA Date"] = pd.to_datetime(tumorVolume["cfDNA Date"])
    tumorVolume["CT Scan Date"] = pd.to_datetime(tumorVolume["CT Scan Date"])
    tumorVolume["cfDNA Days"] = -1000
    tumorVolume["CT Scan Days"] = -1000

    # use C1D1 as day 0
    tumorVolume["cfDNA Days"] = (tumorVolume["cfDNA Date"] - tumorVolume["C1D1"]).dt.days
    tumorVolume["CT Scan Days"] = (tumorVolume["CT Scan Date"] - tumorVolume["C1D1"]).dt.days

    '''
    # use first blood draw date as day 0
    for p in tumorVolume["Patient"].unique():
        baseDate = tumorVolume[(tumorVolume["Patient"] == p) & (tumorVolume["Timepoint"]=="Base")]["cfDNA Date"].values[0]
        tumorVolume.loc[tumorVolume["Patient"] == p, "cfDNA Days"] = (tumorVolume["cfDNA Date"] - baseDate).dt.days
        tumorVolume.loc[tumorVolume["Patient"] == p, "CT Scan Days"] = (tumorVolume["CT Scan Date"] - baseDate).dt.days
    '''
    
    tumorVolume.to_csv("tumor-volume-ctDNA-fraction.tsv", sep="\t", index=False)
    
    

    '''
    # make the slope df

    # both ctDNA and tumor volume must be present
    slope = tumorVolume[(~tumorVolume["ctDNAFraction"].isnull()) & (~tumorVolume["Volume (cm3)"].isnull())].copy()

    slope["cfDNA Date"] = pd.to_datetime(slope["cfDNA Date"])
    slope["CT Scan Date"] = pd.to_datetime(slope["CT Scan Date"])

    slope["ctDNAFraction"] = slope["ctDNAFraction"].diff()
    slope["Volume (cm3)"] = slope["Volume (cm3)"].diff()
    slope["cfDNA Date"] = slope["cfDNA Date"].diff()
    slope["CT Scan Date"] = slope["CT Scan Date"].diff()

    slope = slope.assign(Period=slope.Timepoint.shift(1))
    slope["Period"] = slope["Period"].astype(str) + ":" + slope["Timepoint"].astype(str)
    slope = slope[slope["Timepoint"] != "Base"]

    slope["cfDNA Days"] = slope["cfDNA Date"].dt.days
    slope["CT Scan Days"] = slope["CT Scan Date"].dt.days
    slope["cfDNA slope"] = slope["ctDNAFraction"] / slope["cfDNA Days"]
    slope["CT Scan slope"] = slope["Volume (cm3)"] / slope["CT Scan Days"]


    slope.to_csv("slope.tsv", sep="\t", index=False)


    #sys.exit(0)
    '''


    
    # make plots

    sns.set_theme(style="white")
    patients = tumorVolume.groupby('Patient1')
    for patientName, gr in patients.groups.items():
        
        #if not patientName.startswith("BW-CF-Nivo-003"):
        #   break
        
        print(patientName)
        patient = patients.get_group(patientName).copy()
        
        
        # plot the plot

        fig, ax = plt.subplots(figsize=(10, 10))
    
        plt.figure(dpi=300)
        plt.tight_layout()

        ax = plt.gca()
        
        #rcParams['figure.figsize'] = 10,10

        #if patientName.startswith("BW-CF-Nivo-003"):
        #    rcParams['figure.figsize'] = 20,10


        maxct = patient["ctDNAFraction"].max()
        maxvol = patient["Volume (cm3)"].max()

        plt.ylim([0, maxct*1.1])
        g = sns.lineplot(
            data=patient,
            x="cfDNA Days", y="ctDNAFraction", linewidth=5, marker="o",  markersize=10, color="black"
        )

        g.set_xlabel(xlabel="Days since Baseline blood draw", size=10)
        g.set_ylabel(ylabel="ctDNA Fraction (Black)", size=10)

        ax2 = ax.twinx()
        ax2.set_ylim([0, maxvol*1.1])

        #plt.ylim(0)
        g = sns.lineplot(
            data=patient,
            x="CT Scan Days", y="Volume (cm3)", linewidth=5, marker="o",  markersize=10, color="blue", ax=ax2
        )

       
        
        
        g.set_ylabel(ylabel="Tumor Volume (cm3, Blue)", size=10)
        g.tick_params(labelsize=10)
        g.set_xlabel(xlabel="Days since Baseline blood draw", size=10)
        #axs[1].margins(x=0.05, y=0.05)
        ax.yaxis.grid(False) # Hide the horizontal gridlines
        ax.xaxis.grid(True) # Show the vertical gridlines
    
    
            
        plt.tight_layout()
        plt.savefig(patientName + ".pdf")
        plt.close()
    



















if __name__ == '__main__':
    main()

