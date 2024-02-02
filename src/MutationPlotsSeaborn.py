'''
Created on Jul 15, 2021

@author: patelj1
'''

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import date
from matplotlib import rcParams
from matplotlib.ticker import ScalarFormatter
from pandas.api.types import CategoricalDtype
import sys
import colorcet as cc
import numpy as np



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
    
    #print(unaltered)
    
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
    df["Progression"] = df["Progression"].apply(lambda x: x.replace(" ", "-"))
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
    
def numericTimepoint(tp):
    
    if tp == "Base":
        return int(0)
    
    n = tp[1:]
    return int(n)
    
    
    '''
    #print(tp)
    
    if tp == "Base":
        return int(0)
    elif tp == "W2":
        return int(8)
    elif tp == "W8":
        return int(32)
    elif tp == "W20":
        return int(35)
    elif tp == "W22":
        return int(36)
    elif tp == "W32":
        return int(37)
    elif tp == "W55":
        return int(38)
    elif tp == "W115":
        return int(39)
    elif tp == "W140":
        return int(40)
    else:
        return int(tp[1:])
    '''

    
    '''
    w = int(tp[1:]) - 6
    
    #return int(math.log(w, 1.5))
    return int(w/4) + 6
    '''
    
    '''
    if w <=8:
        return w
    else:
        w = 8 + log2(w - 8)
        return int(w)
    '''



def paddedTimepoint(tp):
    
    if(tp=="Base"):
        return "Base"
    elif(len(tp)==2):
        return tp.replace("W", "W00")
    elif(len(tp)==3):
        return tp.replace("W", "W0")
    else:
        return tp
    







####################################33###########################3


def main():
    
    d = "/Users/patelj1/current/Beryl-17-180/analysis-new"
    
    mutationsFile = d + "/all-with-ctDNAFraction-oncokb-WES.tsv"
    badMutationsFile = d + "/bad-mutations.tsv"
    tumorVolumeFile = d + "/tumor-volume-correlation/tumor-volume-ctDNA-fraction.tsv"
    badMutationsNormalsFile = d + "/../../PromiseStudyJillian/analysis/baseline-8-weeks-ctdna-fraction/bad-mutations/bad-mutations-normals.tsv"
    badMutationsPlasmaFile = d + "/../../PromiseStudyJillian/analysis/baseline-8-weeks-ctdna-fraction/bad-mutations/bad-mutations-plasma.tsv"
    

    
    os.chdir(d + "/seaborn-plots")
    
    
    
    mutations = pd.read_csv(mutationsFile, sep="\t", index_col=False, keep_default_na=False, low_memory=False,
                            dtype={"HOTSPOT": str, "HOTSPOT_INTERNAL": str,
                            "cmo_hotspot": str, "hotspot_whitelist": str, "IS-A-HOTSPOT": str})
    

    # load and prepare tumor volume data
    tumorVolume = pd.read_csv(tumorVolumeFile, sep="\t", index_col=False, keep_default_na=False, low_memory=False)
    tumorVolume = tumorVolume[["Patient", "Volume (cm3)", "CT Scan Days"]]
    tumorVolume = tumorVolume[tumorVolume["CT Scan Days"] != ""]
    #tumorVolume.loc[tumorVolume["CT Scan Days"] == "", "CT Scan Days"] = np.nan
    tumorVolume["CT Scan Days"] = tumorVolume["CT Scan Days"].astype(float)
    tumorVolume["CT Scan Weeks"] = round(tumorVolume["CT Scan Days"]/7)
    tumorVolume["TimepointInt"] = tumorVolume["CT Scan Weeks"]
    tumorVolume["Timepoint"] = tumorVolume["CT Scan Weeks"].apply(lambda x: "Base" if x<=0 else "W" + str(int(x)))
    tumorVolume["Timepoint"] = tumorVolume["Timepoint"].apply(paddedTimepoint)
    tumorVolume = tumorVolume.sort_values(by=["Patient", "Timepoint"])

        
        
    badMutationsNormals = pd.read_csv(badMutationsNormalsFile, sep="\t", index_col=False, keep_default_na=False, low_memory=False)
    badMutationsPlasma = pd.read_csv(badMutationsPlasmaFile, sep="\t", index_col=False, keep_default_na=False, low_memory=False)

    badMutationsNormals["M"] = badMutationsNormals["Mutation"].apply(lambda x: ":".join(x.split(":")[2:6]))
    badMutationsPlasma["M"] = badMutationsPlasma["Mutation"].apply(lambda x: ":".join(x.split(":")[2:6]))
                
    pd.options.display.max_columns = None
    pd.options.display.max_rows = None
    
   
    badMutations = set()
   
    # read bad mutations
    with open(badMutationsFile) as f:
        for line in f:
            words = line.split("\t")
            
            if(words[0]=="Mutation"):
                continue
    
            badMutations.add(words[0])
    
        
    
    
    # remove cfDNA-only mutations
    mutations = mutations[mutations["CalledIn"]!="cfDNA"]
    
    
    mutations["TimepointInt"] = mutations["Timepoint"].apply(numericTimepoint)
    mutations["TimepointInt"] = mutations["TimepointInt"].astype(int)
    mutations["Timepoint"] = mutations["Timepoint"].apply(paddedTimepoint)
    
    #mutations["TimepointNumeric"] = mutations["Timepoint"].apply(numericTimepoint)
    #mutations["TimepointNumeric"] = mutations["TimepointNumeric"].astype(int)
    
    
    # prioritize the mutations for the patient based on:
    # hotspot/oncokb status, called in, synonymous, 
    # sort the df by patient, timepoint and then these factors
    # give color to only first 12? The rest are black with some alpha.
    
    
    mutations["Category"] = "Ordinary"
    
    mutations.loc[(mutations["ONCOGENIC"]=="Likely Oncogenic") | (mutations["ONCOGENIC"]=="Oncogenic")
                  | (mutations["ONCOGENIC"]=="Resistance"), "Category"] = "OncoKB"
                  
    mutations.loc[(mutations["hotspot_whitelist"]=="True") | (mutations["IS-A-HOTSPOT"]=="Y"), "Category"] = "Hotspot"
    
                  
                  
    # remove bad mutations if they are not hotspot
    mutations = mutations[(mutations["Category"]=="Hotspot") | (~mutations["M"].isin(badMutations))]


    # remove additional bad mutations using MSK ACCESS random cohort
    b1 = mutations[(mutations["Category"]!="Hotspot") & (mutations["M"].isin(badMutationsNormals["M"]))]["M"]

    b2 = mutations[(mutations["Category"]!="Hotspot") & (mutations["M"].isin(badMutationsPlasma["M"]))]["M"]

    #print(len(b1))
    #print(len(b2))

    # looks like nothing to remove!!    

    #mutations = mutations[(mutations["Category"]=="Hotspot") | (~mutations["M"].isin(badMutationsNormals["M"]))]

    #mutations = mutations[(mutations["Category"]=="Hotspot") | (~mutations["M"].isin(badMutationsPlasma["M"]))]
    
    
    # remove patient-specific mutations??
    
                  
    mutations["Mutation_Name (CalledIn, Category)"] = mutations["Hugo_Symbol"] + "." + \
                        mutations["HGVSp_Short"] + " (" + mutations["CalledIn"] + ", " + mutations["Category"] + ")"
    
    # remove mutations that have 0 AF in all samples in the same patient
    mutations['MaxAF'] = mutations.groupby('ID')['AF'].transform('max')
    mutations.loc[mutations["MaxAF"]==0, "Mutation_Name (CalledIn, Category)"] = "INVALID"
    
    
    catCategory = CategoricalDtype(
            ['Hotspot', 'OncoKB', 'Ordinary'], ordered=True)
    
             
    catMutationType = CategoricalDtype(
            ['Nonsense_Mutation', 'Frame_Shift_Del', 'Frame_Shift_Ins', 'Missense_Mutation', 
             'Splice_Site', 'In_Frame_Del', 'In_Frame_Ins', 'Silent'], 
                    ordered=True)
    
    catCalledIn = CategoricalDtype(
            ['Both', 'IMPACT', 'cfDNA'], ordered=True)
    
    mutations['Variant_Classification'] = mutations['Variant_Classification'].astype(catMutationType)
    mutations['CalledIn'] = mutations['CalledIn'].astype(catCalledIn)
    mutations['Category'] = mutations['Category'].astype(catCategory)
    
    mutations = mutations.sort_values(by=['Patient', 'Timepoint', 'Category', 
                            'CalledIn', 'Mutation_Name (CalledIn, Category)'])
    
    # save
    mutations.to_csv("all-with-ctDNAFraction-oncokb-processed.tsv" , sep="\t", index=False)


    



    # plot

    '''
    # make the grid plot

    # keep only the patients that have base, w2 and w8
    # 3, 6, 7, 8, 9, 10, 12, 18, 19
    patients = ["BW-CF-Nivo-003", "BW-CF-Nivo-006", "BW-CF-Nivo-007", "BW-CF-Nivo-008", "BW-CF-Nivo-009", "BW-CF-Nivo-010", "BW-CF-Nivo-012", "BW-CF-Nivo-018", "BW-CF-Nivo-019"]

    gmutations = mutations[mutations["Patient"].isin(patients)].copy()
    gmutations = gmutations[gmutations["Timepoint"].isin(["Base", "W002", "W008"])]
    gmutations = gmutations[gmutations['Mutation_Name (CalledIn, Category)'] != "INVALID"]

    sns.set_theme(style="whitegrid")
    plt.tight_layout()

    p = sns.color_palette(cc.glasbey, n_colors=15)

    g = sns.FacetGrid(gmutations, row="Patient", height=5, aspect=2) #hue='Mutation_Name (CalledIn, Category)')
    g.map(sns.lineplot, "Timepoint", "AF", 'Mutation_Name (CalledIn, Category)',palette=p)
    
    for ax in g.axes.ravel():
       ax.legend()
       sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))

    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig("grid.pdf")
    plt.close()
    '''

    

    mutations["ACCESSMSIScore"] = mutations["ACCESSMSIScore"].astype(float)
    tumorVolume["Volume (cm3)"] = tumorVolume["Volume (cm3)"].astype(float)

    zeromuts = []

    # plot the per-patient plots
    sns.set_theme(style="whitegrid")
    patients = mutations.groupby('Patient1')
    for patientName, gr in patients.groups.items():
        
        #if not patientName.startswith("BW-CF-Nivo-003"):
        #    break
        
        print(patientName)
        patient = patients.get_group(patientName)
        
        
        # plot the plots
    
        plt.figure(dpi=300)
        plt.tight_layout()
        
        fig, axs = plt.subplots(3,1,
                          figsize=(20,20),
                          sharex=True,
                          gridspec_kw=dict(height_ratios=[1, 2, 1]))
        
        
        #0
        #palette = {"Breast":"crimson",
        #           "Lung":"green"}
        
        
        maxmsi = patient["ACCESSMSIScore"].max()
        #print(maxmsi)

        g = sns.lineplot(
            data=patient,
            x="TimepointInt", y="ACCESSMSIScore", ci=None, linewidth=5, marker="o",  markersize=20, ax=axs[0]
        )

        axs[0].set_ylim([-0.01, maxmsi*1.1+0.01])
        
        g.set_ylabel(ylabel="ACCESS MSI Score", size=30)
        g.tick_params(labelsize=30)
        g.set(xlabel=None)
        axs[0].margins(x=0.05, y=0.05)
        
        
        
        tp="Base"
        if patientName.startswith("BW-CF-Nivo-011"):
            tp="W26"
        elif patientName.startswith("BW-CF-Nivo-013"):
            tp="W26"
        elif patientName.startswith("BW-CF-Nivo-020"):
            tp="W52"
        elif patientName.startswith("BW-CF-Nivo-021"):
            tp="W19"
        
        patientMutations = patient.loc[patient["Timepoint"]==tp, "Mutation_Name (CalledIn, Category)"]
        patientMutations = patientMutations[patientMutations != "INVALID"]
        numMutations = len(patientMutations)
        
        if numMutations == 0:
            zeromuts.append(patientName)

            g = sns.lineplot(
                data=patient,
                x="TimepointInt", y="ctDNAFraction", linewidth=10, marker="o",  markersize=20, color="black", ax=axs[1]
            )
            
            
            g.set_ylabel(ylabel="AF (ctDNA Fraction in Thick Black)", size=30)
            g.tick_params(labelsize=30)
            g.set_xlabel(xlabel="Timepoint (Weeks since C1D1)", size=30)
            axs[1].margins(x=0.05, y=0.05)
            #axs[1].get_legend().remove()
        
        else:
            
            chosenMutations = patientMutations.tolist()
            #print(len(chosenMutations))
            
            patient = patient[patient["Mutation_Name (CalledIn, Category)"].isin(chosenMutations)]
            
            p = sns.color_palette(cc.glasbey, n_colors=len(chosenMutations))
        
            g = sns.lineplot(
                data=patient,
                x="TimepointInt", y="AF", hue="Mutation_Name (CalledIn, Category)", linewidth=5, marker="o",
                        markersize=20, palette=p, ax=axs[1]
                            )
            
            g = sns.lineplot(
                data=patient,
                x="TimepointInt", y="ctDNAFraction", linewidth=10, marker="o",  markersize=20, color="black", ax=axs[1]
            )
            
            lt = "Mutation_Name (CalledIn, Category) " + "[" + str(numMutations) + "]"
            leg = axs[1].legend(title=lt, fontsize=30, bbox_to_anchor=(1.05, 1))
            
            '''
            handles, labels = axs[1].get_legend_handles_labels()
            newHandles = handles[:20]
            newLabels = labels[:20]
            leg = axs[1].legend(newHandles, newLabels, title=lt, fontsize=30, bbox_to_anchor=(1.05, 1))
            '''
            
            # set the linewidth of each legend object
            # this is so stupid!
            for legobj in leg.legendHandles:
                legobj.set_linewidth(10)
            
            
            # this is also quite stupid!!
            plt.setp(leg.get_title(),fontsize='30')
            
            g.set_ylabel(ylabel="AF (ctDNA Fraction in Thick Black)", size=30)
            g.tick_params(labelsize=30)
            g.set_xlabel(xlabel="Timepoint (Weeks since C1D1)", size=30)
            axs[1].margins(x=0.05, y=0.05)


            p = patientName.split(" ")[0]
            df = tumorVolume[tumorVolume["Patient"]==p]
            if df.empty:
                continue

            maxvolume = df["Volume (cm3)"].max()

            g = sns.lineplot(
                data=df,
                x="TimepointInt", y="Volume (cm3)", linewidth=10, marker="o",  markersize=20, color="blue", ax=axs[2]
            )

            axs[2].set_ylim([-1, maxvolume*1.1])
            
            
            g.set_ylabel(ylabel="Tumor Volume (cm3)", size=30)
            g.tick_params(labelsize=30)
            g.set_xlabel(xlabel="Timepoint (Weeks since C1D1)", size=30)
            axs[2].margins(x=0.05, y=0.05)




        
        
            
        plt.tight_layout()
        plt.savefig(patientName + ".pdf")
        plt.close()


    print(zeromuts)
    

    
    
    


         


if __name__ == '__main__':
    main()
    
    
    
    
    