'''
Created on Jun 10, 2021

@author: patelj1
'''


'''
Make serial timepoint plots and create files for Oncoprinter

'''

import sys
import os 
import pandas as pd
import altair as alt
from datetime import timedelta
from pandas.api.types import CategoricalDtype
import warnings

warnings.simplefilter(action='ignore', category=FutureWarning)





def makeOncoprinterFiles(rowset):
    
    # sort by patient id
    
    df = rowset.copy()
   
    
    
    # fix variant classification
    # Oncoprinter: "MISSENSE", "INFRAME", "TRUNC", "SPLICE", "PROMOTER", or "OTHER" for a mutation alteration
    
    df["Variant_Classification-Original"] = df["Variant_Classification"]
    df["Variant_Classification"] = "OTHER"
    
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
    
    df.loc[df["HGVSp_Short"]=="", "HGVSp_Short"] = "_"
    df["HGVSp_Short"] = df["HGVSp_Short"].apply(lambda x: x.replace(" ", "-"))
    df["BestResponse"] = df["BestResponse"].astype(str).apply(lambda x: x.replace(" ", "-"))
   
    
    # make and write mutation data
    # sort genes by number of occurrences, keep only genes that have been mutated in at least 3 samples
    df['geneFrequency'] = df.groupby('Hugo_Symbol')['Tumor_Sample_Barcode'].transform('nunique')
    #df = df.loc[df["geneFrequency"]>2]
    df = df.sort_values('geneFrequency', ascending=False)
    mutationData = df[["Tumor_Sample_Barcode", "Hugo_Symbol", "HGVSp_Short", "Variant_Classification"]]
    mutationData.to_csv("oncoprinter-mutations.tsv", sep="\t", header=False, index=False)
    
    
    # make and write clinical data
    df = df.drop_duplicates(subset="Tumor_Sample_Barcode")
    
    # this doesn't really determine the oncoprinter sample order
    df = df.sort_values(["Patient", "Timepoint"], 
                                             ascending = (True, True))
    df.rename(columns = {"Tumor_Sample_Barcode": "Sample", "ACCESSMSIScore": "ACCESSMSIScore(number)", 
                         "WESMSIScore": "WESMSIScore(number)"}, inplace=True)
    
    df = df[["Sample", "Timepoint", "BestResponse", "ACCESSMSIScore(number)", "WESMSIScore(number)"]]
    df.to_csv("oncoprinter-clinical.tsv", sep="\t", index=False)
    

    

def main():
    
    d = "/Users/patelj1/current/ACCESSGRAILMetasaticBreastCancer/analysis/"
    
    cfDNAFile = d + "genotypes-access-relevant-annotated.maf"
    tissueFile = d + "genotypes-impact-relevant-annotated.maf"
    manifestFile = d + "manifest-relevant.tsv"
    regimensFile = d + "regimens.tsv"
    responseFile = d + "clinical-response.tsv"
    cfDNACollectionsFile = d + "cfDNA-collection-dates.tsv"
    normalDuplexFile = d + "genotypes-normal-duplex.tsv"
    normalUnfilteredFile = d + "genotypes-normal-unfiltered.tsv"
    impactTumorFile = d + "genotypes-impact-tumor.tsv"
    impactNormalFile = d + "genotypes-impact-normal.tsv"
    badMutationsNormalsFile = d + "bad-mutations-access-cohort/" + "bad-mutations-normals.tsv"
    badMutationsPlasmaFile = d + "bad-mutations-access-cohort/" + "bad-mutations-plasma.tsv"
    badMutationsManualFile = d + "bad-mutations-manual-inspection.tsv"


    
    
    os.chdir(d + "plots-altair-zoomed-in")
    
    
    cfDNAMutations = pd.read_csv(cfDNAFile, sep="\t", index_col=False, keep_default_na=False, low_memory=False)
    tissueMutations = pd.read_csv(tissueFile, sep="\t", index_col=False, keep_default_na=False, low_memory=False)
    manifest = pd.read_csv(manifestFile, sep="\t", index_col=False, keep_default_na=False, low_memory=False)
    normalDuplex = pd.read_csv(normalDuplexFile, sep="\t", index_col=False, keep_default_na=False, low_memory=False)
    normalUnfiltered = pd.read_csv(normalUnfilteredFile, sep="\t", index_col=False, keep_default_na=False, low_memory=False)
    impactTumor = pd.read_csv(impactTumorFile, sep="\t", index_col=False, keep_default_na=False, low_memory=False)
    impactNormal = pd.read_csv(impactNormalFile, sep="\t", index_col=False, keep_default_na=False, low_memory=False)
    badMutationsNormals = pd.read_csv(badMutationsNormalsFile, sep="\t", index_col=False, keep_default_na=False, low_memory=False)
    badMutationsPlasma = pd.read_csv(badMutationsPlasmaFile, sep="\t", index_col=False, keep_default_na=False, low_memory=False)
    badMutationsManual = pd.read_csv(badMutationsManualFile, sep="\t", index_col=False, keep_default_na=False, low_memory=False)
    
    
    # process manifest
    manifest["MRN"] = manifest["MRN"].astype(str)
    manifest = manifest[manifest["Sample_Class"]=="Tumor"]
    manifest["Collection_date"] = pd.to_datetime(manifest["Collection_date"])
    manifest["FirstDate"] = manifest.groupby("Corrected_Investigator_Patient_ID")["Collection_date"].transform("min")
    manifest["Timepoint (Days)"] = (manifest['Collection_date'] - manifest['FirstDate']).dt.days
    #manifest["Timepoint (Weeks)"] = (manifest["Timepoint (Weeks)"]/7).astype(int)
    manifest = manifest[["Corrected_Investigator_Patient_ID", "MRN", "CMO_Sample_ID", "FirstDate", 
                         "Collection_date", "Timepoint (Days)"]]
    manifest = manifest.sort_values(by=["Corrected_Investigator_Patient_ID", "Collection_date"])
    
    manifest.to_csv("manifest.tsv", sep="\t", index=False)
    
    manifest = manifest.drop_duplicates(subset=["Corrected_Investigator_Patient_ID", "Collection_date"])
    
    
    
    # combine cfDNA and tissue mutations
    # common mutations have already been marked CalledIn "Both" in cfDNA file
    # so those can be dropped from tissue file
    cfDNAMutations["ID"] = cfDNAMutations.apply(lambda row: 
                            row["Patient"] + ":" + row["Chromosome"] + ":" + str(row["Start_Position"])
                            + ":" + row["Reference_Allele"] + ":" + row["Tumor_Seq_Allele2"], axis=1)


    cfDNAMutations["CalledIn"] = "cfDNA"
     
    tissueMutations["ID"] = tissueMutations.apply(lambda row: 
                            row["Patient"] + ":" + row["Chromosome"] + ":" + str(row["Start_Position"])
                            + ":" + row["Reference_Allele"] + ":" + row["Tumor_Seq_Allele2"], axis=1)
    
    tissueMutations["CalledIn"] = "Tissue"
    
    
    # remove germline impact mutations
    germline = ["MSK-MB-0057:17:17122440:-:TTCTGTACTCTCTGGCAACACAGGGGCT", "MSK-MB-0042:13:32914438:T:-", 
                "MSK-MB-0037:11:108203626:A:C", "MSK-MB-0073:17:41267744:T:A", "MSK-MB-0044:22:29121230:C:T",
                "MSK-MB-0055:22:29090054:G:A"]
    
    tissueMutations = tissueMutations[~tissueMutations["ID"].isin(germline)]


    # remove TERT promoter mutations other than hotspots
    cfDNAMutations = cfDNAMutations[(cfDNAMutations["Hugo_Symbol"]!="TERT") | ((cfDNAMutations["ID"].str.endswith("5:1295228:G:A")) | (cfDNAMutations["ID"].str.endswith("5:1295250:G:A")))]
    
    # mark common mutations as "Both"
    cfDNAMutations.loc[cfDNAMutations["ID"].isin(tissueMutations["ID"]), "CalledIn"] = "Both"
    
    # also write out common mutations
    common = cfDNAMutations[cfDNAMutations["ID"].isin(tissueMutations["ID"])]
    common.to_csv("called-in-both.maf", sep="\t", index=False)
    
    # add tissue only mutations to cfDNA mutations
    tissueOnly = tissueMutations[~(tissueMutations["ID"].isin(cfDNAMutations["ID"]))]
    #(tissueOnly[["ID"]])
    mutations = pd.concat([cfDNAMutations, tissueOnly], ignore_index=True)
    
    # add normal genotypes
    normalDuplex = normalDuplex.drop("Sample", axis=1)
    normalDuplex = normalDuplex.drop_duplicates(subset=["ID"])   
    mutations = mutations.merge(normalDuplex, how="left", left_on='ID', right_on='ID')
    
    normalUnfiltered = normalUnfiltered.drop("Sample", axis=1)
    normalUnfiltered = normalUnfiltered.drop_duplicates(subset=["ID"])   
    mutations = mutations.merge(normalUnfiltered, how="left", left_on='ID', right_on='ID')

    # add IMPACT genotypes
    impactTumor = impactTumor.drop("Sample", axis=1)
    impactTumor = impactTumor.drop_duplicates(subset=["ID"])   
    mutations = mutations.merge(impactTumor, how="left", left_on='ID', right_on='ID')
    mutations["IMPACT Tumor"] = mutations["IMPACT Tumor"].astype(str)

    impactNormal = impactNormal.drop("Sample", axis=1)
    impactNormal = impactNormal.drop_duplicates(subset=["ID"])   
    mutations = mutations.merge(impactNormal, how="left", left_on='ID', right_on='ID')
    mutations["IMPACT Normal"] = mutations["IMPACT Normal"].astype(str)
    
   
    
    # determine max coverage for each patient-mutation combination
    mutations["maxCoverage"] = mutations.groupby("ID")["Waltz_MD_t_depth"].transform("max")
    
    # remove patient-mutation combinations that have <200 coverage in all samples
    print("Removing patient-mutation combinations that have <200 coverage in all samples of the patient:")
    lowCoverage = mutations[mutations["maxCoverage"]<200]["ID"]
    print(lowCoverage)
    
    mutations = mutations[mutations["maxCoverage"]>=200]



    # remove bad mutations identified at access cohort level
    badMutationsNormals["M"] = badMutationsNormals["Mutation"].apply(lambda x: ".".join(x.split(":")[:2]))
    badMutationsPlasma["M"] = badMutationsPlasma["Mutation"].apply(lambda x: ".".join(x.split(":")[:2]))

    mutations["M"] = mutations["Hugo_Symbol"] + "." + mutations["HGVSp_Short"]

    mutations = mutations[~mutations["M"].isin(badMutationsNormals["M"])]
    mutations = mutations[~mutations["M"].isin(badMutationsPlasma["M"])]

    # remove manually identified bad mutations
    mutations = mutations[~mutations["ID"].isin(badMutationsManual["ID"])]


    # remove mutations in genes that start with "LOC"
    mutations = mutations[~mutations["Hugo_Symbol"].str.startswith("LOC")]
    
    
    
    
    # if read support is 1, make AF 0
    mutations["AF"] = mutations["Waltz_MD_t_alt_count"]/mutations["Waltz_MD_t_depth"]
    mutations.loc[mutations["Waltz_MD_t_alt_count"]==1, "AF"] = 0.0
    
    # make 0 AF values non-zero for the log scale
    mutations.loc[mutations["AF"]==0.0, "AF"] = 0.001

    
    
    # add and rename columns
    mutations["Mutation_Name (CalledIn)"] = mutations["M"] + " (" + mutations["CalledIn"] + ")"
    
    #mutations["Patient1"] = mutations["Patient"] + " (" + mutations["BestResponse"] + ")"
    
    
    # merge mutations and manifest, only keep patients/samples present in the manifest
    mutations = mutations.merge(manifest, how="left", left_on='Tumor_Sample_Barcode', right_on='CMO_Sample_ID')    
    mutations = mutations.dropna(subset=["Collection_date"])
    mutations['Patient (MRN)'] = mutations['Patient'] + " (" + mutations["MRN"] + ")"
       
    mutations = mutations.sort_values(by=["Patient", "Timepoint (Days)", "Mutation_Name (CalledIn)"])
    

    
    # identify CHIP mutations
    chip = mutations[(mutations["Hugo_Symbol"]!="RB1") & (mutations["Hugo_Symbol"]!="PTEN")].copy()
    chip["Duplex Normal"] = chip["Duplex Normal"].astype(str)
    chip["Unfiltered Normal"] = chip["Unfiltered Normal"].astype(str)

    
    chip["ds"] = chip["Duplex Normal"].apply(lambda x: x.split(" ")[0].split("/")[0])
    chip["dt"] = chip["Duplex Normal"].apply(lambda x: x.split(" ")[0].split("/")[1])
    chip["ds"] = chip["ds"].astype(int)
    chip["dt"] = chip["dt"].astype(float)
    chip["df"] = chip["ds"]/chip["dt"]

    chip["us"] = chip["Unfiltered Normal"].apply(lambda x: x.split(" ")[0].split("/")[0])
    chip["ut"] = chip["Unfiltered Normal"].apply(lambda x: x.split(" ")[0].split("/")[1])
    chip["us"] = chip["us"].astype(int)
    chip["ut"] = chip["ut"].astype(float)
    chip["uf"] = chip["us"]/chip["ut"]

    #chip = mutations.copy()

    chip = chip[((chip["ds"]>=2) & (chip["df"]>=0.005)) | ((chip["us"]>=3) & (chip["uf"]>=0.005))]

    chip.to_csv("chip.tsv", sep="\t", index=False)

    mutations = mutations[~mutations["ID"].isin(chip["ID"])]



    # remove patient-mutation combinations where CalledIn is cfDNA and maxAF is 0.001
    mutations["maxAF"] = mutations.groupby("ID")["AF"].transform("max")
    print("Removing patient-mutation combinations where CalledIn is cfDNA and maxAF is 0.001:")
    zeroAFcfDNAOnly = mutations[(mutations["CalledIn"]=="cfDNA") & (mutations["maxAF"]<=0.001)]["ID"]
    print(zeroAFcfDNAOnly)
    
    mutations = mutations[~mutations["ID"].isin(zeroAFcfDNAOnly)]


    # process oncokb and info
    mutations["IS-A-HOTSPOT"] = mutations["IS-A-HOTSPOT"].astype(str)
    mutations["ONCOGENIC"] = mutations["ONCOGENIC"].astype(str)
    mutations["OncoKB"] = mutations["ONCOGENIC"]
    mutations.loc[mutations["IS-A-HOTSPOT"]=="Y", "OncoKB"] = mutations["OncoKB"] + ", HOTSPOT "

    
    '''
    # plot only RB1, PTEN, ESR1, ERBB2, PIK3CA, BRCA1, BRCA2
    #interestingGenes = ["BRCA2"]
    interestingGenes = ["RB1", "PTEN", "ESR1", "ERBB2", "PIK3CA", "BRCA1", "BRCA2"]
    
    mutations = mutations[mutations["Hugo_Symbol"].isin(interestingGenes)]
    os.chdir("../plots-altair-interesting")
    '''

    # save
    mutations.to_csv("all.tsv", sep="\t", index=False)

    
    
    # load regimens
    regimens = {}
    startDates = {}
    endDates = {}
    
    f = open(regimensFile)
    
    count = 0
    for line in f:
        count += 1
        
        if count == 1:
            continue
        
        words = line.strip().split("\t")
        p = words[0]
        r = words[1].split(",")
        s = words[3].split(",")
        e = words[4].split(",")
        
        r1 = []
        s1 = []
        e1 = []
        for i, item in enumerate(r):
            if item == "+" or item == "":
                continue
            
            # arrange regimen combination lexicographically
            x = r[i].split("+")
            x.sort()
            r[i] = "+".join(filter(None, x))
            
            
            r1.append(r[i])
            s1.append(s[i])
            e1.append(e[i])
        
        regimens[p] = r1
        startDates[p] = s1
        endDates[p] = e1
        
    #print(regimens.keys())
    f.close()
    
    
    # load clinical response
    response = pd.read_csv(responseFile, sep="\t", index_col=False, keep_default_na=False, low_memory=False)
    response["Date"] = pd.to_datetime(response["Date"])
    
    # add fixed colors
    #1f77b4, #ff7f0e, #2ca02c, #d62728, #9467bd,
    #8c564b, #e377c2, #7f7f7f, #bcbd22, #17becf
    
    responseValues = ["Complete Response", "Mixed Response", "Partial Response", "Progression of Disease (Mixed Response)",
                        "Progression of Disease", "Stable Disease", "Stable Disease/Partial Response"]
    
    responseColors = ["#2ca02c", "#9467bd", "#1f77b4", "#e377c2", "#d62728", "#ff7f0e", "#8c564b" ]
    
    
    
    # load cfDNA collections dates
    cfDNACollections = pd.read_csv(cfDNACollectionsFile, sep="\t", index_col=False, keep_default_na=False, low_memory=False)
    cfDNACollections["cfDNA Collection Date"] = pd.to_datetime(cfDNACollections["cfDNA Collection Date"])
    
    # add fixed colors
    collectionValues = ["Sequenced", "Not Sequenced"]
    collectionColors = ["#2ca02c", "#d62728"]
    

    
    # make the plots

    patients = mutations.groupby('Patient')
    for patientName, samples in patients.groups.items():
        
        #if not patientName.startswith("MSK-MB-0007"):
        #    continue
        
        
        r = regimens[patientName]
        s = startDates[patientName]
        e = endDates[patientName]
        
        # make dataframe for shading
        rdf = pd.DataFrame({
            'Regimen': r,
            'StartDate': s,
            'EndDate': e
            })
        
        # convert dates to timepoints
        rdf['StartDate'] = pd.to_datetime(rdf['StartDate'])
        rdf['EndDate'] = pd.to_datetime(rdf['EndDate'])
        f = manifest[manifest["Corrected_Investigator_Patient_ID"]==patientName]['FirstDate'].iloc[0]
        rdf["StartTimepoint"] = (rdf['StartDate'] - f).dt.days
        #regimens["StartTimepoint"] = (regimens["StartTimepoint"]/7).astype(int)
        rdf["EndTimepoint"] = (rdf['EndDate'] - f).dt.days
        #regimens["EndTimepoint"] = (regimens["EndTimepoint"]/7).astype(int)
    
       
        
        
        print(patientName)
        patient = patients.get_group(patientName)
        patientNameMRN = patient["Patient (MRN)"].iloc[0]
        
        patientResponse = response[response["Patient"]==patientName]
        patientcfDNACollections = cfDNACollections[cfDNACollections["Patient"]==patientName]

        #if len(patientcfDNACollections.index) == 0:
        #    continue

        minDate = patientcfDNACollections["cfDNA Collection Date"].min() - timedelta(days=30)
        maxDate = patientcfDNACollections["cfDNA Collection Date"].max() + timedelta(days=30)
        
        
        chart = alt.Chart(
            rdf.reset_index()
        ).mark_rect(
            opacity=0.1
        ).encode(
            x = alt.X('StartDate', scale=alt.Scale(domain=[minDate, maxDate])),
            #x='StartDate',
            x2='EndDate',
            y=alt.value(0),  # pixels from top
            y2=alt.value(300),  # pixels from top
            #color='Regimen'
            color=alt.Color('Regimen', sort=['StartDate']),
            tooltip=["Regimen", "StartDate", "EndDate"]
        ).properties(title = patientNameMRN, width=400, height=300).interactive()
        
        
        
        chart += alt.Chart(patient).mark_line().encode(
            x=alt.X('Collection_date', title="Time"),
            #x=alt.X('Timepoint (Days)'),
            y=alt.Y('AF', scale=alt.Scale(type='log', base=10)),
            color=alt.Color('Mutation_Name (CalledIn)', legend=None),
            
        ) 
        
        chart += alt.Chart(patient).mark_point(size=50).encode(
            x=alt.X('Collection_date', title=""),
            #x=alt.X('Timepoint (Days)'),
            y=alt.Y('AF'), # scale=alt.Scale(type='log', base=10)),
            color=alt.Color('Mutation_Name (CalledIn)'),
            #shape='cfDNAExclusive',
            tooltip=["Patient", "Collection_date", "ID", "Mutation_Name (CalledIn)", "OncoKB", "AlphaMissense", "AF", "Waltz_MD_t_depth", 
                      "Waltz_MD_t_alt_count", "Duplex Normal", "Unfiltered Normal", "IMPACT Tumor", "IMPACT Normal"],
           
        ).interactive()
        
        chart = chart.resolve_scale(color='independent')
        
        
        
        chart1 = alt.Chart(patientResponse).mark_circle(size=15).encode(
             x=alt.X('Date', title="Clinical Response"),
            #x=alt.X('Timepoint (Days)'),
            #y=alt.Y('AF'), #scale=alt.Scale(type='log', base=10)),
            color=alt.Color('Clinical Response', scale=alt.Scale(domain=responseValues, range=responseColors)),
            #shape='cfDNAExclusive',
            tooltip=["Patient", "Date", "Clinical Response"],
        ).interactive()
        
        chart1 = chart1.resolve_scale(color='independent')
        
        
        
        chart2 = alt.Chart(patientcfDNACollections).mark_circle(size=15).encode(
            x=alt.X('cfDNA Collection Date', title="cfDNA Collection"), 
            color=alt.Color('cfDNA Collection Status', scale=alt.Scale(domain=collectionValues, range=collectionColors)),
            tooltip=["Patient", "cfDNA Collection Date", "cfDNA Collection Status"],
        ).interactive()
        
        chart2 = chart2.resolve_scale(color='independent')
        
                
        combined = alt.vconcat(chart, chart1, chart2).resolve_scale(x="shared")
        
        combined.save(patientName + ".html")
    
    
    
    
    
    #######
    



if __name__ == '__main__':
    main()
    
    
    
    
    