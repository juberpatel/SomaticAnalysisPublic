'''
Author: Juber Patel

'''

import os
import pandas as pd
from scipy import stats
import sys


def computeGeometricMean(s):
    
    nonZero = s[s!=0]
    
    if nonZero.count() == 0:
        return 0
    else:
        return stats.gmean(nonZero)

'''
Compute ctDNA fractions using cfDNA data only, ie without using tissue data
'''
def main():
    
    d = "/Users/patelj1/current/PromiseStudyJillian/analysis/baseline-8-weeks-ctdna-fraction"

    os.chdir(d)

    cfDNAFile = "merged-AFGM-oncokb.tsv"

    patientsFile = "patients-with-baseline-8-weeks-jillian.tsv"

    diseaseFile = "disease-treatment-arm.tsv"

    badPlasmaMutationsFile = "bad-mutations/bad-mutations-plasma.tsv"

    badNormalsMutationsFile = "bad-mutations/bad-mutations-normals.tsv"

    manualRemovalFile = "manual-removal.tsv"

    goodTERTMutationsFile = "impact-frequent-TERT-mutations.tsv"



    cfDNA = pd.read_csv(cfDNAFile, sep="\t", index_col=False, keep_default_na=False, low_memory=False)


    # make 1 read support 0 read support
    cfDNA.loc[cfDNA["AD"] == 1, "AD"] = 0

    # calculate AF
    cfDNA["AF"] = cfDNA["AD"]/cfDNA["DP"]

    cfDNA["Mutation"] = cfDNA.apply(lambda row: ":".join([row["Hugo_Symbol"], row["HGVSp_Short"], row["M"]]), axis=1)
    cfDNA["PatientMutation"] = cfDNA.apply(lambda row: ":".join([row["Patient"], row["Mutation"]]), axis=1)

    # restrict to patients of interest
    patients = pd.read_csv(patientsFile, sep="\t", index_col=False, keep_default_na=False, low_memory=False)

    cfDNA = cfDNA[cfDNA["Patient"].isin(patients["Patient"])]

    print(len(cfDNA.index))

    # remove bad mutations
    badPlasmaMutations = pd.read_csv(badPlasmaMutationsFile, sep="\t", index_col=False, keep_default_na=False, low_memory=False)

    badPlasmaMutations["M"] = badPlasmaMutations["Mutation"].apply(lambda x: ":".join(x.split(":")[2:6]))    

    badNormalsMutations = pd.read_csv(badNormalsMutationsFile, sep="\t", index_col=False, keep_default_na=False, low_memory=False)

    badNormalsMutations["M"] = badNormalsMutations["Mutation"].apply(lambda x: ":".join(x.split(":")[2:6]))

    cfDNA = cfDNA[~cfDNA["M"].isin(badPlasmaMutations["M"])]
    cfDNA = cfDNA[~cfDNA["M"].isin(badNormalsMutations["M"])]

    # remove genes starting with "LOC"
    cfDNA = cfDNA[~cfDNA["Hugo_Symbol"].str.startswith("LOC")]

    # remove manual removals
    manualRemovals = pd.read_csv(manualRemovalFile, sep="\t", index_col=False, keep_default_na=False, low_memory=False)

    cfDNA = cfDNA[~cfDNA["PatientMutation"].isin(manualRemovals["Manual"])]

    # only keep good TERT mutations
    goodTERTMutations = pd.read_csv(goodTERTMutationsFile, sep="\t", index_col=False, keep_default_na=False, low_memory=False)

    badTERTMutations = cfDNA[(cfDNA["Hugo_Symbol"]=="TERT") & (~cfDNA["M"].isin(goodTERTMutations["M"]))]

    cfDNA = cfDNA[~cfDNA["M"].isin(badTERTMutations["M"])]
    
    print(len(cfDNA.index))

    


    # adjust timepoints
    cfDNA = cfDNA[cfDNA["Timepoint"] != "Progression"]
    cfDNA.loc[cfDNA["Timepoint"] == "8 Weeks-Progression", "Timepoint"] = "8 Weeks"

    print(len(cfDNA.index))

    # if a mutation has 0 AF in both timepoints, remove it
    cfDNA["MaxAF"] = cfDNA.groupby(["PatientMutation"])["AF"].transform("max")
    cfDNA = cfDNA[cfDNA["MaxAF"] > 0.0]
    # drop max AF column
    cfDNA = cfDNA.drop(columns=["MaxAF"])

    
    # only keep mutations called in at least 1 timepoint
    cfDNA["CalledInPatient"] = cfDNA.groupby(["PatientMutation"])["Called"].transform(lambda x: "Yes" in x.values)
    cfDNA = cfDNA[cfDNA["CalledInPatient"] == True]

    # mark hotspot mutations
    cfDNA["Hotspot"] = cfDNA.apply(lambda row: row["IS-A-HOTSPOT"]=="Y" or row["ONCOGENIC"]=="Oncogenic", axis=1)


    '''
    approach 0:
    # compute median in baseline, remove mutations with AF < median
    t = cfDNA.copy()
    t = t[t["Timepoint"] == "Baseline"]
    t["MedianAF"] = t.groupby(["Patient"])["AF"].transform("median")
    t = t[t["AF"] >= t["MedianAF"]]
    cfDNA = cfDNA[cfDNA["PatientMutation"].isin(t["PatientMutation"])]
    '''

    '''
    # approach 1: use baseline timepoint
    # determine clonal mutations at baseline
    t = cfDNA.copy()
    t = t[t["Timepoint"] == "Baseline"]
    t["HighestAF"] = t.groupby(["Patient"])["AF"].transform("max")
    t = t[t["AF"] >= t["HighestAF"] * 0.3]
    cfDNA = cfDNA[cfDNA["PatientMutation"].isin(t["PatientMutation"])]

        
    # compute median AF
    cfDNA["MedianAF"] = cfDNA.groupby(["Patient", "Timepoint"])["AF"].transform("median")
    '''

    '''
    # approach 2:
    if the patient has hotspot mutations
        pick the timepoint where hotspot mutations have higher AF??
    else
        Pick the timepoint where mutations have higher AF??
    '''

    hotspots = cfDNA.copy()
    hotspots = hotspots[hotspots["Hotspot"] == True]
    hotspots["SelectedTimepoint"] = "Undecided"

    
    #hotspots["MedianAF"] = hotspots.groupby(["Patient", "Timepoint"])["AF"].transform("mean")
    #hotspots["MedianAF"] = hotspots.groupby(["Patient", "Timepoint"])["AF"].transform("median")
    hotspots["MedianAF"] = hotspots.groupby(["Patient", "Timepoint"])["AF"].transform(lambda x: computeGeometricMean(x))

    # use the timepoint with higher median for identifying clonal mutations
    hotspots = hotspots.sort_values(by=["Patient", "MedianAF"], ascending=[True, False])
    hotspots = hotspots.drop_duplicates(subset=["Patient"], keep="first")
    hotspots["SelectedTimepoint"] = hotspots["Timepoint"]

    cfDNA = cfDNA.merge(hotspots[["Patient", "SelectedTimepoint"]], on=["Patient"], how="left")

    #cfDNA.to_csv("t.tsv", sep="\t", index=False)
    #sys.exit(0)

    '''
    if there are crosstrends AND trendsum.abs() < 0.01 (difficult cases)
        then try to decide based on oncogenic mutations
    
    if this doesnt break the tie, just pick baseline my man
    '''

     # if you see crosstrends in hotspot mutations, give preference to oncogenic mutations
    hotspots = cfDNA.copy()
    hotspots = hotspots.drop(columns=["SelectedTimepoint"])
    hotspots = hotspots[hotspots["Hotspot"] == True]
    hotspots["NumMutations"] = hotspots.groupby(["Patient", "Timepoint"])["AF"].transform("count")
    hotspots = hotspots[hotspots["NumMutations"] > 1]

    # find trend for each patient-mutation
    # b -a = b + (-a)
    hotspots["T"] = hotspots["AF"]
    hotspots.loc[hotspots["Timepoint"] == "Baseline", "T"] = -hotspots["T"]
    hotspots["Trend"] = hotspots.groupby(["PatientMutation"])["T"].transform('sum')
    hotspots["Trendsum"] = hotspots.groupby(["Patient"])["Trend"].transform('sum')
    hotspots["Trendsum"] = hotspots["Trendsum"]/2

    # remove cases with good trendsum
    hotspots = hotspots[hotspots["Trendsum"].abs() < 0.01]

    # remove patients where trend is always positive or always negative
    hotspots["PositiveTrend"] = hotspots.groupby(["Patient"])["Trend"].transform(lambda x: (x > 0).all())
    hotspots["NegativeTrend"] = hotspots.groupby(["Patient"])["Trend"].transform(lambda x: (x < 0).all())
    hotspots = hotspots[(hotspots["PositiveTrend"] == False) & (hotspots["NegativeTrend"] == False)]

    hotspots["HasOncogenic"] = hotspots.groupby(["Patient"])["ONCOGENIC"].transform(lambda x: "Oncogenic" in x.values)

    # if a patient has oncogenic mutations, then pick the timepoint where oncogenic mutations have higher AF
    t = hotspots.copy()
    t = t[t["ONCOGENIC"] == "Oncogenic"]
    t["OncogenicAF"] = t.groupby(["Patient", "Timepoint"])["AF"].transform("max")
    t = t.sort_values(by=["Patient", "OncogenicAF"], ascending=[True, False])
    t = t.drop_duplicates(subset=["Patient"], keep="first")
    t["SelectedTimepoint"] = t["Timepoint"]
    hotspots = hotspots.merge(t[["Patient", "SelectedTimepoint"]], on=["Patient"], how="left")

    # if a patient has no oncogenic mutations, then pick baseline as selected timepoint  
    hotspots.loc[hotspots["HasOncogenic"] == False, "SelectedTimepoint"] = "Baseline"
    hotspots = hotspots.drop_duplicates(subset=["Patient"], keep="first")
    hotspots["SelectedTimepointDifficult"] = hotspots["SelectedTimepoint"]

    cfDNA = cfDNA.merge(hotspots[["Patient", "SelectedTimepointDifficult"]], on=["Patient"], how="left")
    cfDNA.loc[~cfDNA["SelectedTimepointDifficult"].isnull(), "SelectedTimepoint"] = cfDNA["SelectedTimepointDifficult"]
    

    rest = cfDNA.copy()
    rest = rest[rest["SelectedTimepoint"].isnull()]

    print("Patients with no hotspot mutations: ")
    print(rest["Patient"].unique())

    #rest["MedianAF"] = rest.groupby(["Patient", "Timepoint"])["AF"].transform("mean")
    #rest["MedianAF"] = rest.groupby(["Patient", "Timepoint"])["AF"].transform("median")
    rest["MedianAF"] = rest.groupby(["Patient", "Timepoint"])["AF"].transform(lambda x: computeGeometricMean(x))
    rest = rest.sort_values(by=["Patient", "MedianAF"], ascending=[True, False])
    rest = rest.drop_duplicates(subset=["Patient"], keep="first")
    rest["SelectedTimepointRest"] = rest["Timepoint"]

    cfDNA = cfDNA.merge(rest[["Patient", "SelectedTimepointRest"]], on=["Patient"], how="left")

    cfDNA.loc[cfDNA["SelectedTimepoint"].isnull(), "SelectedTimepoint"] = cfDNA["SelectedTimepointRest"]

    # determine clonal mutations in the selected timepoint
    t = cfDNA.copy()
    t = t[t["Timepoint"] == t["SelectedTimepoint"]]

    # determine clonal mutations
    t["HighestAF"] = t.groupby(["Patient"])["AF"].transform("max")
    t = t[t["AF"] >= t["HighestAF"] * 0.3]
    cfDNA = cfDNA[cfDNA["PatientMutation"].isin(t["PatientMutation"])]

    # compute median AF
    #cfDNA["MedianAF"] = cfDNA.groupby(["Patient", "Timepoint"])["AF"].transform("mean")
    #cfDNA["MedianAF"] = cfDNA.groupby(["Patient", "Timepoint"])["AF"].transform("median")
    cfDNA["MedianAF"] = cfDNA.groupby(["Patient", "Timepoint"])["AF"].transform(lambda x: computeGeometricMean(x))

    cfDNA["MedianAF"] = cfDNA["MedianAF"] * 2
    
    cfDNA.to_csv("cleaned-cfDNA.tsv", sep="\t", index=False)

    print(len(cfDNA["Patient"].unique()))

    df = cfDNA.copy()
    df = df[['Patient', 'MedianAF', 'Timepoint']]
    df = df.drop_duplicates(subset=["Patient", "Timepoint"], keep="first")

    # add missing patients
    missingPatients = patients[~patients["Patient"].isin(df["Patient"])]["Patient"].unique()

    print(missingPatients)

    l = []
    for p in missingPatients:
        l.append([p, 0.0, "Baseline"])
        l.append([p, 0.0, "8 Weeks"])

    missing = pd.DataFrame(l, columns=df.columns)

    df = pd.concat([df, missing], ignore_index=True)
    
    
    # add disease info
    disease = pd.read_csv(diseaseFile, sep="\t", index_col=False, keep_default_na=False, low_memory=False)

    disease = disease[["ID", "Disease", "Treatment Arm"]]
    disease.columns = ["Patient", "Disease", "TreatmentArm"]

    df = df.merge(disease, on="Patient", how="left")


    df.to_csv("cfDNA-fractions.tsv", sep="\t", index=False)
    















if __name__ == '__main__':
    main()
