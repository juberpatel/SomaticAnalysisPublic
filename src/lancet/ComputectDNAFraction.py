import os
import pandas as pd
from scipy import stats




'''
use

t_alt = (v[TCN*P + NCN*(1−P)])/P

check if this value is closer to mcn or lcn 
'''
def determineMutantCN(row):
    
    lcn = row["lcn"]
    
    '''
    if lcn == "NA" or lcn == "0":
        lcn = "1"
    '''

    if lcn == "NA":
        lcn = "0"

    
    lcn = int(lcn)
    tcn = int(row["tcn"])
    mcn = tcn - lcn
    
    # if mcn and lcn are equal, we are done
    if mcn == lcn:
        row["mutant_cn"] = lcn

    # if lcn is 0, we are done
    elif lcn == 0:
        row["mutant_cn"] = mcn
    
    else:
        ncn = 2

        if row["Chromosome"] == "X" and row["Sex"] == "M":
            ncn = 1

        purity = row["Purity"]
        
        vaf = row["Tumor_MAF"]
        t_alt = vaf * (tcn * purity + ncn * (1.0 - purity))
        t_alt = t_alt/purity
        row["computed_t_alt_cn"] = t_alt
        
        dist1 = abs(t_alt - mcn)
        dist2 = abs(t_alt - lcn)
        
        if dist1 < dist2:
            row["mutant_cn"] = mcn
        else:
            row["mutant_cn"] = lcn
    
        #print(row["Patient"] + "\t" + str(mcn) + "\t" + str(lcn) + "\t" + str(vaf) + "\t" + 
        #str(purity) + "\t" + str(t_alt) + 
                                                    #  "\t" + str(row["mutant_cn"]))  
    
    return row



'''

VAF = (TALT * P)/(TCN * P + NCN * (1−P))

Therefore,

P = (NCN * VAF)/(TALT + (NCN−TCN) * VAF)

Purity in cfDNA sample is the ctDNA fraction


where

VAF = variant allele fraction (expected)
P = tumor purity
TALT = alternate copies in tumor
TCN = total copies in tumor
NCN = total copies in normal (2)

'''
def computeCTDNAFraction(row):

    #print(row["PatientMutation"])
    
    ncn = 2

    if row["Chromosome"] == "X" and row["Sex"] == "M":
            ncn = 1

    vaf = row["AF"]
    t_alt = row["mutant_cn"]
    tcn = int(row["tcn"])
    
    #print(row["Tumor_Sample_Barcode"] + "\t" + row["Mutation"] + "\t" + str(ncn) + "\t" + str(vaf) + "\t" + str(t_alt) + "\t" + str(tcn))
    
    ctDNAFraction = (ncn * vaf)/(t_alt + ((ncn-tcn) * vaf))
    
    row["MutationctDNAFraction"] = ctDNAFraction
    
    return row

########################################################################   

def computeGeometricMean(s):
    
    nonZero = s[s!=0]
    
    if nonZero.count() == 0:
        return 0
    else:
        return stats.gmean(nonZero)



#####################################################################


'''

compute ctDNA fraction for each cfDNA sample

1. Use only clonal tissue mutations
2. Compute the mutant allele copy number for each clonal mutation
3. Compute ctDNA fraction for each clonal mutation
4. For each cfDNA sample, take geometric mean of the ctDNA fractions of clonal mutations in that sample
5. The tisue mutations file should have only filtered clonal mutations and must have these fields: 
            Patient, mutation, TUMOR_MAF, tcn, lcn, purity


'''
def main():
    
    d = "/Users/patelj1/current/PromiseStudyJillian/analysis/baseline-8-weeks-ctdna-fraction"

    os.chdir(d)

    #impactFile = "baseline-8-weeks-ABSOLUTE.tsv"
    impactFile = "baseline-8-weeks-clonal-ABSOLUTE.tsv"

    absoluteFile = "ABSOLUTE-summary.tsv"
    genotypingFile = "genotypes.maf"
    timepointsFile = "timepoints.tsv"
    diseaseFile = "disease-treatment-arm.tsv"



    # process impact data
    impact = pd.read_csv(impactFile, sep="\t", index_col=False, keep_default_na=False, low_memory=False)

    # turn 23 to X
    impact["Chromosome"] = impact["Chromosome"].astype(str)
    impact.loc[impact["Chromosome"] == "23", "Chromosome"] = "X"

    # remove P093 from consideration, we don't have proper tissue sample
    impact = impact[impact["Patient"] != "P093"]

    # some bad samples were rescued with manual inspection
    # remove bad samples
    badSamples = ["P-0032151-T02-IM6", "P-0046437-T05-IM7", "P-0047891-T03-IM7"]

    impact = impact[~impact["Tumor_Sample_Barcode"].isin(badSamples)]

    impact["Mutation"] = impact.apply(lambda x: x["Hugo_Symbol"] + x["HGVSp_Short"] + ":" + str(x["Chromosome"]) + ":" + str(x["Start_Position"]) + ":" + x["Reference_Allele"] + ":" + x["Tumor_Seq_Allele2"], axis=1)
    impact["Mutation"] = impact["Mutation"].str.replace("p.", ".")

    impact["PatientMutation"] = impact.apply(lambda x: x["Patient"] + ":" + x["Mutation"], axis=1)

    impact["tcn"] = impact["tcn"].astype(str)
    impact["lcn"] = impact["lcn"].astype(str)


    absolute = pd.read_csv(absoluteFile, sep="\t", index_col=False, keep_default_na=False, low_memory=False)
    absolute = absolute[["sample", "purity", "ploidy"]]
    absolute.columns = ["Tumor_Sample_Barcode", "Purity", "Ploidy"]

    impact = pd.merge(impact, absolute, on="Tumor_Sample_Barcode", how="left")

    # set certain tumor purity values manually
    impact.loc[impact["Tumor_Sample_Barcode"] == "P-0045924-T01-IM6", "Purity"] = 0.1
    impact.loc[impact["Tumor_Sample_Barcode"] == "P-0035487-T01-IM6", "Purity"] = 0.1
    impact.loc[impact["Tumor_Sample_Barcode"] == "P-0050711-T01-IM6", "Purity"] = 0.1
    # P044
    impact.loc[impact["Tumor_Sample_Barcode"] == "P-0000001-T01-IM6", "Purity"] = 0.6
    

    # use clonal mutations from all samples
    # but remove duplicate mutations based on purity

    impact = impact.sort_values(by=["PatientMutation", "Purity"], ascending=[True, False])
    impact = impact.drop_duplicates(subset=["PatientMutation"], keep="first")

    '''
    # from each patient, keep the sample with highest purity
    impact = impact.sort_values(by=["Patient", "Purity", "Tumor_Sample_Barcode"], ascending=[True, False, True])
    samplesToKeep = impact.drop_duplicates(subset=["Patient"], keep="first")["Tumor_Sample_Barcode"]
    impact = impact[impact["Tumor_Sample_Barcode"].isin(samplesToKeep)]
    '''

    impact["Tumor_MAF"] = impact.apply(lambda x: float(x["t_alt_count"])/(float(x["t_ref_count"]) + float(x["t_alt_count"])), axis=1)

    impact.to_csv("tissue.tsv", sep="\t", index=False)

    

    # process genotyping data
    genotypes = pd.read_csv(genotypingFile, sep="\t", index_col=False, keep_default_na=False, low_memory=False)

    genotypes["Tumor_Sample_Barcode"] = genotypes["Tumor_Sample_Barcode"].apply(lambda x: x.split("_")[0])
    genotypes["Patient"] = genotypes["Tumor_Sample_Barcode"].apply(lambda x: x.split("-")[0])

    genotypes["Mutation"] = genotypes.apply(lambda x: x["Hugo_Symbol"] + x["HGVSp_Short"] + ":" + str(x["Chromosome"]) + ":" + str(x["Start_Position"]) + ":" + x["Reference_Allele"] + ":" + x["Tumor_Seq_Allele2"], axis=1)
    genotypes["Mutation"] = genotypes["Mutation"].str.replace("p.", ".")

    genotypes["PatientMutation"] = genotypes.apply(lambda x: x["Patient"] + ":" + x["Mutation"], axis=1)

    # only keep patient-mutation combinations that are clonal in impact data
    genotypes = genotypes[genotypes["PatientMutation"].isin(impact["PatientMutation"])]
    
    # remove normal samples
    genotypes["SampleType"] = genotypes["Tumor_Sample_Barcode"].apply(lambda x: "Normal" if x.split("-")[1].startswith("N") else "Tumor")
    genotypes = genotypes[genotypes["SampleType"] == "Tumor"]

     # attach timepoint info
    timepoints = pd.read_csv(timepointsFile, sep="\t", index_col=False, keep_default_na=False, low_memory=False)
    timepoints.columns = ["Tumor_Sample_Barcode", "Timepoint"]
    genotypes = pd.merge(genotypes, timepoints, on="Tumor_Sample_Barcode", how="left")

    # remove Progression timepoint and harmonize 8-week timepoint
    genotypes = genotypes[genotypes["Timepoint"] != "Progression"]
    genotypes.loc[genotypes["Timepoint"]=="8 Weeks-Progression", "Timepoint"] = "8 Weeks"

    p1 = genotypes["Patient"].unique()
    print(len(p1))

    # remove mutations not covered in ACCESS
    genotypes["minCoverage"] = genotypes.groupby("PatientMutation")["Waltz_MD_t_depth"].transform("min")
    genotypes = genotypes[genotypes["minCoverage"] >= 200]

    p2 = genotypes["Patient"].unique()
    print(len(p2))

    print("Not covered in ACCESS:")
    print(set(p1) - set(p2))


    # calculate AF
    genotypes.loc[genotypes["Waltz_MD_t_alt_count"] == 1, "Waltz_MD_t_alt_count"] = 0
    genotypes["AF"] = genotypes["Waltz_MD_t_alt_count"] / genotypes["Waltz_MD_t_depth"]

    # remove mutations with AF == 0 across timepoints
    #genotypes["maxAF"] = genotypes.groupby("PatientMutation")["AF"].transform("max")
    #genotypes = genotypes[genotypes["maxAF"] != 0]

    

    # attach disease and treatment arm info
    disease = pd.read_csv(diseaseFile, sep="\t", index_col=False, keep_default_na=False, low_memory=False)
    disease = disease[["ID", "Disease", "Treatment Arm", "Sex"]]
    disease.columns = ["Patient", "Disease", "TreatmentArm", "Sex"]
    genotypes = pd.merge(genotypes, disease, on="Patient", how="left")

    # attach impact data
    impact = impact[["PatientMutation", "Tumor_MAF", "Purity", "Ploidy", "tcn", "lcn"]]
    genotypes = pd.merge(genotypes, impact, on="PatientMutation", how="left")

    genotypes.to_csv("cfDNA.tsv", sep="\t", index=False)



    
    # compute ctDNA fraction

    chosenMutations = genotypes.copy()

    # determine mutant allele copy number
    chosenMutations = chosenMutations.apply(determineMutantCN, axis=1)
    
    
    # compute ctDNA fraction   
    chosenMutations = chosenMutations.apply(computeCTDNAFraction, axis=1)
    
    # get geometric mean as ctDNAFraction
    chosenMutations["ctDNAFraction"] = chosenMutations["MutationctDNAFraction"].groupby(
                    chosenMutations["Tumor_Sample_Barcode"]).transform(computeGeometricMean)
    
    chosenMutations.to_csv("t3.tsv", sep="\t", index=False)
    
    # get the useful info involed in calculation of ctDNA fraction
    chosenMutations = chosenMutations[["Tumor_Sample_Barcode", "Mutation", "computed_t_alt_cn", "mutant_cn", "MutationctDNAFraction", "ctDNAFraction"]]
    
    genotypes = genotypes.merge(chosenMutations, how="left", on=["Tumor_Sample_Barcode", "Mutation"])
    
    genotypes.to_csv("baseline-8-weeks-ctDNA-fractions.tsv", sep="\t", index=False)

    genotypes = genotypes[["Patient", "Timepoint", "Disease", "TreatmentArm", "ctDNAFraction"]]
    genotypes = genotypes.drop_duplicates(subset=["Patient", "Timepoint"], keep="first")
    
    genotypes.to_csv("ctDNA-fractions.tsv", sep="\t", index=False)
    
    







if __name__ == "__main__":
    main()