'''

Author: Juber Patel

'''


import os
import pandas as pd
from scipy import stats


def main():
    
    d = "/Users/patelj1/current/PromiseStudyJillian/analysis/baseline-8-weeks-ctdna-fraction/bad-mutations"

    os.chdir(d)

    randomPlasmaFile = "genotypes-random-plasma.maf"
    randomNormalsFile = "genotypes-random-normals.maf"


    df = pd.read_csv(randomPlasmaFile, sep="\t", index_col=False, keep_default_na=False, low_memory=False)

    df["Mutation"] = df.apply(lambda row: row["Hugo_Symbol"] + ":" + row["HGVSp_Short"] + ":" + row["Chromosome"] + ":" + str(row["Start_Position"]) + ":" + row["Reference_Allele"] + ":" + row["Tumor_Seq_Allele2"], axis=1)

    df = df[df["Waltz_total_t_depth"]>=50]
    df.loc[df["Waltz_total_t_alt_count"]==1, "Waltz_total_t_alt_count"] = 0
    df["AF"] = df["Waltz_total_t_alt_count"]/df["Waltz_total_t_depth"]
    
    
    df["CoveredIn"] = df.groupby("Mutation")["AF"].transform("count")
    df["PresentIn"] = df.groupby("Mutation")["AF"].transform(lambda x: len(x[x>0]))
    df["PresentInFraction"] = df["PresentIn"]/df["CoveredIn"]
    df["MeanAF"] = df.groupby("Mutation")["AF"].transform("mean")
    df["MedianAF"] = df.groupby("Mutation")["AF"].transform("median")
    df["MinAF"] = df.groupby("Mutation")["AF"].transform("min")
    df["MaxAF"] = df.groupby("Mutation")["AF"].transform("max")


    df = df[["Mutation", "CoveredIn", "PresentIn", "PresentInFraction", "MeanAF", "MedianAF", "MinAF", "MaxAF"]]
    df = df.drop_duplicates(subset=["Mutation"])
    df = df.sort_values(by=["PresentInFraction"], ascending=False)
    df = df.round(5)

    df.to_csv("presence-in-random-plasma.tsv", sep="\t", index=False)


    df = pd.read_csv(randomNormalsFile, sep="\t", index_col=False, keep_default_na=False, low_memory=False)

    df["Mutation"] = df.apply(lambda row: row["Hugo_Symbol"] + ":" + row["HGVSp_Short"] + ":" + row["Chromosome"] + ":" + str(row["Start_Position"]) + ":" + row["Reference_Allele"] + ":" + row["Tumor_Seq_Allele2"], axis=1)

    df = df[df["Waltz_total_t_depth"]>=50]
    df.loc[df["Waltz_total_t_alt_count"]==1, "Waltz_total_t_alt_count"] = 0
    df["AF"] = df["Waltz_total_t_alt_count"]/df["Waltz_total_t_depth"]
    
    
    df["CoveredIn"] = df.groupby("Mutation")["AF"].transform("count")
    df["PresentIn"] = df.groupby("Mutation")["AF"].transform(lambda x: len(x[x>0.005]))
    df["PresentInFraction"] = df["PresentIn"]/df["CoveredIn"]
    df["MeanAF"] = df.groupby("Mutation")["AF"].transform("mean")
    df["MedianAF"] = df.groupby("Mutation")["AF"].transform("median")
    df["MinAF"] = df.groupby("Mutation")["AF"].transform("min")
    df["MaxAF"] = df.groupby("Mutation")["AF"].transform("max")


    df = df[["Mutation", "CoveredIn", "PresentIn", "PresentInFraction", "MeanAF", "MedianAF", "MinAF", "MaxAF"]]
    df = df.drop_duplicates(subset=["Mutation"])
    df = df.sort_values(by=["PresentInFraction"], ascending=False)
    df = df.round(5)

    df.to_csv("presence-in-random-normals.tsv", sep="\t", index=False)










if __name__ == '__main__':
    main()
