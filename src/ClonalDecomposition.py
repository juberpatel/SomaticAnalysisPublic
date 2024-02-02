'''
Created on Oct 27, 2020

@author: Juber Patel
'''

import pandas as pd
import os
import seaborn as sb
import copy
import matplotlib.pyplot as plt
import statistics
from collections import OrderedDict
import sys




def makeMYIDOld(gene, protein, ID, hotspot1, hotspot2, hotspot3, hotspot4, oncogenic, highestLevel):
    myid = str(gene) + "." + str(protein) + ";"
    
    if oncogenic=="Oncogenic":
        myid = myid + "Oncogenic"
    elif oncogenic=="Likely Oncogenic":
        myid = myid + "LikelyOncogenic"
    elif oncogenic=="Predicted Oncogenic":
        myid = myid + "PredictedOncogenic"
    
    
    if highestLevel!="":
        words = highestLevel.split("LEVEL_")
        myid = myid + ":" + words[1]
        
    myid = myid + ";"
    
    if hotspot1==True or hotspot2==True or hotspot3==True or hotspot4=="Y":
        myid = myid + "HOTSPOT"
        
    myid = myid + ";" + ID
    
    return myid



def makeMutationBlocksFile(patientName, m, df):
    
    identifyEWPSamples = open(patientName + "-mutation-blocks-for-signatures.txt", "w")
    identifyEWPSamples.write("Sample\tConstituents\tMutation\tChr\tPosition\tRef\tAlt\n")
    
    groups = m.groupby(["presentIn"]).groups
    for group, labels in groups.items():
        indices = [i for i, val in enumerate(group) if val=="1"]
        samples = [m.columns[i] for i in indices]
        for label in labels:
            line = "block-" + group + "\t" + ",".join(samples)
            record = df[df["MYID"]==label].iloc[0]
            line = line + "\t" + record["MYID"] + "\t" + str(record["Chromosome"]) + \
                        "\t" + str(record["Start_Position"]) + "\t" \
                    + record["Reference_Allele"] + "\t" + record["Tumor_Seq_Allele2"]
                    
            identifyEWPSamples.write(line + "\n")



def makeBranchPrivateMutationsFile(patientName, m, df):
    
    # write trunk mutations.
    # per sample, write whole sample, branch and private mutations
    
    f = open(patientName + "-branch-private-mutations-for-signatures.txt", "w")
    f.write("Sample\tMutation\tChr\tPosition\tRef\tAlt\n")
    
    numSamples = len(m.iloc[0]["presentIn"])
    groups = m.groupby(["presentIn"]).groups
    
    #trunk
    trunk = []
    trunkGroup = "1"*numSamples
    for label in groups[trunkGroup]:
        record = df[df["MYID"]==label].iloc[0]
        line = record["MYID"] + "\t" + str(record["Chromosome"]) + "\t" + str(record["Start_Position"]) + "\t" \
                + record["Reference_Allele"] + "\t" + record["Tumor_Seq_Allele2"] + "\n"
        trunk.append(line)
        
    # write trunk mutations
    for line in trunk:
        f.write("trunk\t" + line)
        
        
    # now do sample by sample
    for i in range(numSamples):
        sample = m.columns[i]
        # collect private mutations
        private = []
        privateGroup = "0"*numSamples
        privateGroup = privateGroup[:i] + "1" + privateGroup[i+1:]
        if privateGroup in groups:
            for label in groups[privateGroup]:
                record = df[df["MYID"]==label].iloc[0]
                line = record["MYID"] + "\t" + str(record["Chromosome"]) + "\t" + str(record["Start_Position"]) + "\t" \
                        + record["Reference_Allele"] + "\t" + record["Tumor_Seq_Allele2"] + "\n"
                private.append(line)
            
        # collect branch mutations ie mutations also present in other samples but not in trunk
        branch = []
        for group, labels in groups.items():
            if group == trunkGroup or group == privateGroup or group[i] != "1":
                continue
            
            for label in labels:
                record = df[df["MYID"]==label].iloc[0]
                line = record["MYID"] + "\t" + str(record["Chromosome"]) + "\t" + str(record["Start_Position"]) + "\t" \
                        + record["Reference_Allele"] + "\t" + record["Tumor_Seq_Allele2"] + "\n"
                branch.append(line)
                
        # write the whole sample
        for line in trunk:
            f.write(sample + "-all\t" + line)
        for line in private:
            f.write(sample + "-all\t" + line)
        for line in branch:
            f.write(sample + "-all\t" + line)
        
        # write private mutations
        for line in private:
            f.write(sample + "-private\t" + line)
            
        # write branch mutations
        for line in branch:
            f.write(sample + "-branch\t" + line)
        
        
def makePhylogicNDTInputFiles(patientName, m, df, absoluteReviewedFolder):
    
    if not os.path.isdir("PhylogicNDT-input"):
        os.mkdir("PhylogicNDT-input")
    
    # read purity file
    summary = pd.read_csv(absoluteReviewedFolder + "/absolute.test.ABSOLUTE.table.txt", 
                               sep="\t", index_col=False, low_memory=False, dtype=str)
    # write sif file
    sif = open("PhylogicNDT-input/" + patientName + ".sif", "w")
    sif.write("sample_id\tmaf_fn\tseg_fn\tpurity\ttimepoint\n")
    
    # make the header 
    header = "Hugo_Symbol\tChromosome\tStart_position\tReference_Allele\tTumor_Seq_Allele2\tt_ref_count\tt_alt_count"
    
    for i in range(0, 100):
        if i < 10:
            header += ("\tccf_0.0" + str(i))
        else:
            header += ("\tccf_0." + str(i))
    
    header += "\tccf_1.00"
    
    print("Generating PhylogicNDT input files for:")

        
    timePoint = 0
    for sample in m.columns:
        if sample == "numSamples" or sample == "sums" or sample == "presentIn":
            continue
        
        print(sample)
        purity = summary[summary["sample"]==sample].iloc[0]["purity"]
        path = "PhylogicNDT-input/" + sample + "_ABS_MAF.txt"
        sif.write(sample + "\t" + "/mnt/" + os.path.basename(os.getcwd()) + "/" + path + 
                  "\t\t" + purity + "\t" + str(timePoint) + "\n")
        timePoint += 1
        
        
        absolute = pd.read_csv(absoluteReviewedFolder + "/SEG_MAF/" + sample + "_ABS_MAF.txt", 
                               sep="\t", index_col=False, low_memory=False, dtype=str)
        
        PNDTInputFile = open(path, "w")
        PNDTInputFile.write(header + "\n")
        
        for mutation in m.index:
            
            dfRecord = df[(df["MYID"]==mutation) & (df["Tumor_Sample_Barcode"]==sample)]
            
            # no df record for this combination of sample and mutation
            # still create a record with 0 CCF (ccf_0=1)
            if(len(dfRecord.index)==0):
                dfRecord = df[(df["MYID"]==mutation)].iloc[0]
                line = dfRecord["Hugo_Symbol"] + "_" + dfRecord["HGVSp_Short"] + \
                        "\t" + str(dfRecord["Chromosome"]) + "\t" + str(dfRecord["Start_Position"]) + \
                        "\t" + dfRecord["Reference_Allele"] + "\t" + dfRecord["Tumor_Seq_Allele2"] + \
                    "\t" + "500" + "\t" + "0"
                
                line += "\t1"
                for i in range(0, 100):
                    line += "\t0"
                
                PNDTInputFile.write(line + "\n")
                continue
            
            dfRecord = dfRecord.iloc[0]
            
            absoluteRecord = absolute[(absolute["sample"]==dfRecord["Tumor_Sample_Barcode"]) &
                                      (absolute["Hugo_Symbol"]==(dfRecord["Hugo_Symbol"] + "_" + dfRecord["HGVSp_Short"])) &
                                      (absolute["Chromosome"]==str(dfRecord["Chromosome"])) &
                                      (absolute["Start_position"]==str(dfRecord["Start_Position"]))]
            
            # make the output line
            line = dfRecord["Hugo_Symbol"] + "_" + dfRecord["HGVSp_Short"] + "\t" + str(dfRecord["Chromosome"]) + "\t" \
                    + str(dfRecord["Start_Position"]) + "\t" + dfRecord["Reference_Allele"] + "\t" + \
                    dfRecord["Tumor_Seq_Allele2"] + "\t" + str(dfRecord["t_ref_count"]) + "\t" + \
                    str(dfRecord["t_alt_count"])
                    
            # if no absolute record, create a record with 0 CCF (ccf_0=1)
            if(len(absoluteRecord.index)==0):
                line += "\t1"
                for i in range(0, 100):
                    line += "\t0"
                
            else:
                absoluteRecord = absoluteRecord.iloc[0].astype("str")
            
                for i in range(0, 100):
                    if i==0:
                        line += ("\t" + absoluteRecord["CCF_0"])
                    elif i < 10:
                        line += ("\t" + absoluteRecord["CCF_0.0" + str(i)])
                    elif i%10==0:
                        line += ("\t" + absoluteRecord["CCF_0." + str(i//10)])
                    else:
                        line += ("\t" + absoluteRecord["CCF_0." + str(i)])
                
                line += ("\t" + absoluteRecord["CCF_1"])
            
            
            
            PNDTInputFile.write(line + "\n")
        
            
        PNDTInputFile.close()
    
    sif.close()



##########################################################################3

def getMYID(mutations, row):
    
    mutation = str(row["Chromosome"]) + ":" + str(row["Start_Position"]) + ":" + \
                        row["Reference_Allele"] + ":" + row["Tumor_Seq_Allele2"]
    
    if mutation in mutations:
        return mutations[mutation]
    else:
        return "----"
  

def makeMYID(mutations, geneNames, row):
    
    '''
    gene.protein;oncokb;hotspot;genomicchange;rs;cosm
    '''
    
    mutation = str(row["Chromosome"]) + ":" + str(row["Start_Position"]) + ":" + \
                        row["Reference_Allele"] + ":" + row["Tumor_Seq_Allele2"]
    
    if mutation in mutations:
        return mutations[mutation]
    
    
    # make the id
    myid = str(row["Hugo_Symbol"]) + "." + str(row["HGVSp_Short"]) + ";"
    
    oncogenic = row["ONCOGENIC"]
    if oncogenic=="Oncogenic":
        myid = myid + "Oncogenic"
    elif oncogenic=="Likely Oncogenic":
        myid = myid + "LikelyOncogenic"
    elif oncogenic=="Predicted Oncogenic":
        myid = myid + "PredictedOncogenic"
    
    
    if row["HIGHEST_LEVEL"]!="":
        words = row["HIGHEST_LEVEL"].split("LEVEL_")
        myid = myid + ":" + words[1]
        
    myid = myid + ";"
    
    if row["HOTSPOT"]==True or row["HOTSPOT_INTERNAL"]==True or row["cmo_hotspot"]==True or row["IS-A-HOTSPOT"]=="Y":
        myid = myid + "HOTSPOT"
        
    myid = myid + ";"
    
    myid = myid + mutation + ";"
    
    rs = ""
    cosm = ""
    words = row["ID"].split(";")
    for word in words:
        if(rs=="" and word.startswith("rs")):
            rs = word
        elif(cosm=="" and word.startswith("COSM")):
            cosm = word
            
    myid = myid + rs + ";"
    myid = myid + cosm
    
    
    
    # adjust mutation labels by cancer geneNames and by mutation class
    gene = row["Hugo_Symbol"]
    if gene in geneNames and row["Variant_Classification"] != "Silent" and row["Variant_Classification"] != "Splice_Region":
        myid = myid.replace(gene, geneNames[gene]) + ";Cancer"
    else:
        myid = myid + ";"
        
    
    mutations[mutation] = myid
    
    return myid


def createPairtreeInputs(patientName, df, tr):
    
    print("Creating Pairtree input...")
    
    # convert tr to dict of total reads
    trDict = {}
    for index, row in tr.iterrows():
        sample = row["Tumor_Sample_Barcode"]
        mutation = row["MYID"]
        dp = row["TUMOR_DP"]
        
        if sample not in trDict:
            trDict[sample] = {}
            
        trDict[sample][mutation] = dp
    
    
    df = df[df["Cluster"]!=-1]
    df = df.sort_values(by=["Cluster"])
    
    
    if not os.path.isdir("Pairtree-input"):
        os.mkdir("Pairtree-input")
    
    # write ssm file
    ssm = open("Pairtree-input/" + patientName + ".ssm", "w")
    ssm.write("id\tname\tvar_reads\ttotal_reads\tvar_read_prob\n")
    
    
    mutationClusters = OrderedDict()
    count = 0
    samples = df["Tumor_Sample_Barcode"].unique()
    # sample index
    sampleIndex = {}
    for i in range(len(samples)):
        sampleIndex[samples[i]] = i
        
    mutationGroups = df.groupby("MYID", sort=False)
    for mutation, group in mutationGroups:
        
        # get total reads from trDict
        total_reads = [0] * len(samples)
        for i in range(len(samples)):
            if mutation in trDict[samples[i]]:
                total_reads[i] = int(trDict[samples[i]][mutation])
            else:
                total_reads[i] = 100
        
        var_reads = [0] * len(samples)
        var_read_prob = [0.5] * len(samples)
        
        cluster = -1
        for index, row in group.iterrows():
            i = sampleIndex[row["Tumor_Sample_Barcode"]]
            var_reads[i] = round(row["TUMOR_MAF"] * row["TUMOR_DP"])
            var_read_prob[i] = float(row["mutant_cn"])/float(row["tcn"])
            # adjust probability for sample purity
            var_read_prob[i] = var_read_prob[i]  * float(row["purity"])
            cluster = row["Cluster"]
            
        
        id = "s" + str(count)
        total_reads = ','.join(map(str, total_reads))
        var_reads = ','.join(map(str, var_reads))
        var_read_prob = ','.join(map(str, var_read_prob))
        
        ssm.write(id + "\t" + mutation + "\t" + var_reads + "\t" + total_reads + "\t" + var_read_prob + "\n")
        
        if cluster not in mutationClusters:
            mutationClusters[cluster] = []
        
        mutationClusters[cluster].append(id)
            
        
        count += 1
        
     
    ssm.close()
    
    # build the params.json line
    #line = '{"samples": ["'
    #line = line + '", "'.join(samples) + '"]'
    
    line = '{"samples": ' + str(samples.tolist()) + ', "clusters": ['
    for cluster in mutationClusters:
        l = mutationClusters[cluster]
        
        print(cluster, "\t", len(l))
        line = line + str(l) + ', '
    
    line = line[:-2] + '], "garbage": []}'
    line = line.replace("'", '"')
    
    # write params.json file
    json = open("Pairtree-input/" + patientName + ".params.json", "w")
    json.write(line + "\n")
    json.close()
    
    print("Wrote " + str(len(mutationClusters)) + " clusters.")
    
    
'''
use

t_alt = (v[TCN*P + NCN*(1−P)])/P

check if this value is closer to mcn or lcn 
'''
def determineMutantCN(lcn, tcn, purity, vaf):

    lcn = int(lcn)
    tcn = int(tcn)
    mcn = tcn - lcn
    
    # if mcn and lcn are equal, we are done
    if mcn == lcn:
        computed_t_alt_cn = lcn
        mutant_cn = lcn
    
    else:
        ncn = 2
        purity = float(purity)
        
        t_alt = vaf * (tcn * purity + ncn * (1.0 - purity))
        t_alt = t_alt/purity
        computed_t_alt_cn = t_alt
        
        dist1 = abs(t_alt - mcn)
        dist2 = abs(t_alt - lcn)
        
        if dist1 < dist2:
            mutant_cn = mcn
        else:
            mutant_cn = lcn
    
    # it can't be 0. Assume the copy number was not called properly and make it 1.    
    if mutant_cn == 0:
        mutant_cn = 1
    
    return computed_t_alt_cn, mutant_cn



def addPurityCN(row, purity, totalCN, minorCN):
    
    row["purity"] = purity[row["Tumor_Sample_Barcode"]]
    purity = row["purity"]
    vaf = row["TUMOR_MAF"]
    
    sample = row["Tumor_Sample_Barcode"]
    chr = str(row["Chromosome"])
    start = row["Start_Position"]
    
    segs = totalCN[sample][chr]
    
    tcn = "2"
    lcn = "1"
    for key in segs:
        s = key.split("-")
        segStart = int(s[0])
        segEnd = int(s[1])
        
        if start >= segStart and start <= segEnd:
            tcn = segs[key]
            lcn = minorCN[sample][chr][key]
            if pd.isna(lcn):
                lcn = "1"
            break
    
    # the rare case where tcn is 0, we assume it is 1 since we are detecting a mutation there.
    if tcn=="0":
        tcn="1"
        
    row["tcn"] = tcn
    row["lcn"] = lcn
    
        
    # add mutant allele copy number
    # just getting the values since determineMutantCN messes up column order!!!
    a, b = determineMutantCN(lcn, tcn, purity, vaf)
    row["computed_t_alt_cn"] = a
    row["mutant_cn"] = b
    
    return row
    
    
   

def loadCNDict(copyNumbers, totalCN, minorCN):
    
    for index, row in copyNumbers.iterrows():
        sample = row["ID"].split("_")[0]
        chr = row["chrom"]
        seg = row["loc.start"] + "-" + row["loc.end"]
        
        if sample not in totalCN:
            totalCN[sample] = {}
            
        if chr not in totalCN[sample]:
            totalCN[sample][chr] = {}
            
        totalCN[sample][chr][seg] = row["tcn.em"]
        
        if sample not in minorCN:
            minorCN[sample] = {}
            
        if chr not in minorCN[sample]:
            minorCN[sample][chr] = {}
            
        minorCN[sample][chr][seg] = row["lcn.em"]
    




###########################################################################33

    
'''
1. Do the needed filtering
2. Make CCF heatmap text file and the plot
3. Assign mutations to clusters
4. Compute mutant allele copy number
5. Generate Pairtree input files

'''       
    
def main():
    
    patient = "MSK-AB-0001"
    d = "/Users/patelj1/current/Autopsy/analysis/filtered-new/" + patient
    os.chdir(d)
    mutationsFile = patient + "-ABSOLUTE-oncokb.tsv"
    totalReadsFile = d + "/../" + patient + ".genotyped_mutations.jlab_and_cmo.filtered.maf-exac-annoated.tsv"
    patientName = patient + "-absolute"
    badSamplesFile = patient + "-bad-samples.txt"
    copyNumberFile = patient + "-copy-numbers.txt"
    absoluteSummaryFile = d + "/ABSOLUTE/Output/absolute_results/Output/reviewed/absolute.test.ABSOLUTE.table.txt"
    
    #samplesThreshold = 12
    mutationsThreshold = 10
    

    
    geneListFile = "/Users/patelj1/current/workspace/SomaticAnalysis/cancer-HRD-impact505-genes.txt"
    mutationTypesFile = "/Users/patelj1/current/workspace/SomaticAnalysis/mutation-types.txt"
    sampleLabelsFile = "/Users/patelj1/current/Autopsy/analysis/labels.txt"
    
    
    # ************************************************
    
    df = pd.read_csv(mutationsFile, sep="\t", index_col=False, keep_default_na=False, low_memory=False,
                     dtype={"HOTSPOT": bool, "HOTSPOT_INTERNAL": bool,
                            "cmo_hotspot": bool, "HIGHEST_LEVEL": str, "Hugo_Symbol": str, "HGVSp_Short": str,
                            "Variant_Classification": str})
    
    # load total reads values from the genotyping file 
    tr = pd.read_csv(totalReadsFile, sep="\t", index_col=False, keep_default_na=False, low_memory=False,
                     dtype={"HOTSPOT": bool, "HOTSPOT_INTERNAL": bool,
                            "cmo_hotspot": bool, "HIGHEST_LEVEL": str, "Hugo_Symbol": str, "HGVSp_Short": str,
                            "Variant_Classification": str})
    tr.loc[tr["Chromosome"]=="X", "Chromosome"] = "23"
    

    # remove mutations present in mouse normal, if MOUSE column is present
    if "MOUSE" in df.columns:
        df = df[df["MOUSE"]!="PRESENT"]

    # this can be changed to another column if necessary
    ccfColumn = "Cancer_Cell_Fraction"
    #ccfColumn = "ccf"   
    df = df[df[ccfColumn]!="NA"]
    # convert ccf column to float after replacing "." values
    df[ccfColumn] = df[ccfColumn].replace(".", "0").astype("float64")
    
    # load cancer gene names 
    geneNames = {}
    with open(geneListFile) as f:
        for line in f:
            gene = line.split()[0]
            a = gene.split("-HRD")[0]
            b = a.split("-IMPACT")[0]
            geneNames[b] = gene
            
    # create a mapping between a genomic mutation ie chr, pos, ref, alt and it's MYID so that the same mutation always 
    # has the same MYID
    mutationIds = {}
    df["MYID"] = df.apply(lambda row: makeMYID(mutationIds, geneNames, row), axis=1)  
    tr["MYID"] = tr.apply(lambda row: getMYID(mutationIds, row), axis=1)
    
    # drop the specified bad samples  
    badSamples = []
    
    with open(badSamplesFile) as f:
        for line in f:
            badSamples.append(line.split()[0])
    
    print("Removing Samples:")
    print(badSamples)
    df = df[~(df["Tumor_Sample_Barcode"].isin(badSamples))]
    #print(df["Tumor_Sample_Barcode"].unique())
    
    
    # remove mutations that are not "called" in at least one remaining sample
    # Jlab.called_or_genotyped and CMO.called_or_genotyped
    df["called1"] = df.groupby('MYID')['Jlab.called_or_genotyped'].transform(lambda x: ",".join(x))
    df["called2"] = df.groupby('MYID')['CMO.called_or_genotyped'].transform(lambda x: ",".join(x))

    
    n1 = len(df.index)
    
    df = df.loc[(df["called1"].str.contains("called")) | (df["called2"].str.contains("called"))]
    
    n2= len(df.index)
    
    print("Mutations before: " + str(n1) + ", Mutations after: " + str(n2))
    
    df = df.drop(['called1', 'called2'], axis=1)
   
    
    
    # read absolute summary file file
    absoluteSummary = pd.read_csv(absoluteSummaryFile, 
                               sep="\t", index_col=False, low_memory=False, dtype=str)
    purity = dict(zip(absoluteSummary["sample"], absoluteSummary["purity"]))
    
    
    # read copy number file
    copyNumbers = pd.read_csv(copyNumberFile, 
                               sep="\t", index_col=False, low_memory=False, dtype=str)
    totalCN = {}
    minorCN = {}
    
    loadCNDict(copyNumbers, totalCN, minorCN)
    
    
    # add purity, tcn, lcn and mutant allele copy number
    df = df.apply(addPurityCN, args=(purity, totalCN, minorCN), axis=1)
    
    
    # use provided sample labels
    sampleLabels = {}
    with open(sampleLabelsFile) as f:
        for line in f:
            tokens = line.split()
            sampleLabels[tokens[0]] = tokens[1]
            
    #print(sampleLabels)
    
    df = df.replace({"Tumor_Sample_Barcode": sampleLabels})
    tr = tr.replace({"Tumor_Sample_Barcode": sampleLabels})

    # only keep needed columns from tr
    tr = tr[["Tumor_Sample_Barcode", "MYID", "TUMOR_DP"]]
   
    # make the CCF matrix 
    ccfMatrix = df.pivot_table(index="MYID", columns="Tumor_Sample_Barcode", values=ccfColumn, fill_value=0)
    
    
    # drop any row (mutation) where all values are 0
    numSamples = ccfMatrix.astype(bool).sum(axis=1)
    ccfMatrix["numSamples"] = numSamples
    ccfMatrix = ccfMatrix[ccfMatrix["numSamples"]>0]
    ccfMatrix.drop(columns="numSamples", inplace=True)
    
    
    # sort on number of samples a mutation is present in, then the exact samples it is present in,
    # then the sum of CCF values in those samples
    numSamples = ccfMatrix.astype(bool).sum(axis=1)
    presentIn = ccfMatrix.astype(bool).replace({True:"1", False:"0"}).apply("".join, axis=1)
    rank = presentIn.apply(lambda s: int(s, 2))
    sums = ccfMatrix.sum(axis=1)

    # add these columns AFTER calculating them outside of the dataframe
    ccfMatrix["numSamples"] = numSamples
    ccfMatrix["sums"] = sums
    ccfMatrix["presentIn"] = presentIn
    ccfMatrix["presentIn"] = ccfMatrix["presentIn"].astype(str)
    #ccfMatrix["rank"] = rank
    #ccfMatrix["rank"] = ccfMatrix["rank"].astype(int)
    
    # do the sorting
    #ccfMatrix.sort_values(by=["rank", "sums"], ascending=False, inplace=True)
    ccfMatrix.sort_values(by=["numSamples", "presentIn", "sums"], ascending=False, inplace=True)
    
    
    # drop a mutation if it is present in only one sample (private) and CCF is < 0.5
    #ccfMatrix = ccfMatrix.loc[(ccfMatrix["numSamples"]>1) | (ccfMatrix["sums"]>=0.5)]
    
    
    # make mutation blocks files for signature analysis
    #makeMutationBlocksFile(patientName, ccfMatrix, df)
    #makeBranchPrivateMutationsFile(patientName, ccfMatrix, df)
    #makePhylogicNDTInputFiles(patientName, ccfMatrix, df, absoluteReviewedFolder)
    
    
    # attach variant classification column to ccfMatrix
    vc = dict(zip(df.MYID, df.Variant_Classification))
    vcList = []
    for myid in ccfMatrix.index.values:
        vcList.append(vc[myid])
        
    ccfMatrix["Variant_Classification"] = vcList
        
        
    
    
    
    # remove mutation blocks (mutations present in exact same samples)
    # if there are 3 or less mutations in the block and none of them is annotated (as above)
    print("Removing Mutation Blocks:")
    groups = ccfMatrix.groupby(["presentIn"], sort=False).groups
    for group, labels in groups.items():
        # values used: 1: 30, 2,3,6,4: 10
        if len(labels)<mutationsThreshold:  #and group.count('1')<samplesThreshold:  
            print(group, "\t", labels)
            ccfMatrix = ccfMatrix[ccfMatrix["presentIn"]!=group]
    
    
    # get cluster number for each mutation
    count = 1
    mutationClusters = {}
    groups = ccfMatrix.groupby(["presentIn"], sort=False).groups
    for group, labels in groups.items():
        print(group + "\t" + str(len(labels)))
        for label in labels:
            mutationClusters[label] = count
            
        count += 1
    
    # add cluster number to df
    df["Cluster"] = df.apply(lambda row: mutationClusters.get(row["MYID"], -1), axis=1)
    
    
    # save matrix and mutations to files
    ccfMatrix = ccfMatrix.drop(["numSamples", "sums", "presentIn"], axis=1)
    t = ccfMatrix.drop(["Variant_Classification"], axis=1).T
    absoluteFile = patientName + "-filtered-annotated-ccfs.txt"
    t.to_csv(absoluteFile, sep="\t")
    
    clusteredFile = mutationsFile.replace(".tsv", "-clustered.tsv")
    df.to_csv(clusteredFile, sep="\t", index=False)
    
    
    # this call changes df!
    createPairtreeInputs(patientName.replace("-absolute", ""), df, tr)
    
    
    '''
    # create dict of possible mutation types
    mutationTypes = {}
    with open(mutationTypesFile) as f:
        next(f)
        for line in f:
            words = line.split()
            mutationTypes[words[0]] = float(words[1])
    
    # create the mutation types matrix using ccfMatrix
    mtMatrix = ccfMatrix.copy()
    vcCol = mtMatrix.columns.get_loc("Variant_Classification")
    
    # replace mutation type strings with code
    for row in range(0, len(mtMatrix.index)):
        value = mtMatrix.iat[row, vcCol]
        if value == "":
            value = 1.0
        else:
            value = mutationTypes[value]
        
        # if ccf was 0, mark it as absent
        for column in range(0, len(mtMatrix.columns)-1):
            if mtMatrix.iat[row, column]==0.0:
                mtMatrix.iat[row, column] = 0.9
            else:
                mtMatrix.iat[row, column] = value - 0.1
    
    
    # save the mutation types matrix
    mtMatrix = mtMatrix.drop(["Variant_Classification"], axis=1).T   
    mtMatrix.to_csv(patientName + "-filtered-annotated-mutation-types.txt", sep="\t")
    '''
    
    
    print("Plotting ABSOLUTE CCF heatmap...")
    absolute = pd.read_csv(absoluteFile, sep="\t", index_col=0, low_memory=False)
    #absolute = absolute.rename(columns = lambda x: x if x=="" else x.split(";")[3])
    
    # reverse column order
    #absolute = absolute[absolute.columns[::-1]]
    
    
    # get cluster numbers
    clusterNumbers = []
    for mut in absolute.columns:
        clusterNumbers.append(mutationClusters[mut])
        
    print(clusterNumbers)
    
    
    fig = plt.figure(dpi=300)
    plt.tight_layout()
    
    
    c = sb.color_palette("Reds", as_cmap=True)
    c = copy.copy(c)
    c.set_under("white")
    
    # make cluster marker colors
    #columnPalette = sb.color_palette(["black", "white"])
    columnPalette = sb.color_palette()
    columnColors = {}
    for mut in mutationClusters:
        indicator = mutationClusters[mut] % 2
        columnColors[mut] = columnPalette[indicator]
        
    columnColors = pd.Series(columnColors)
    
    r = sb.clustermap(absolute, method="complete", metric="canberra", col_cluster=False, 
                      cmap=c, vmin=0.01, vmax=1.0, cbar_kws={'label': 'CCF'}, figsize=(20, 12),
                      xticklabels=True, yticklabels=True, col_colors=[columnColors])
    r.ax_row_dendrogram.set_visible(False)
    #plt.setp(r.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    
    labels = []
    for tick in r.ax_heatmap.xaxis.get_major_ticks():
        label = tick.label1
        label.set_size(15)
        text = label.get_text()
        name = text.split(";")[0]
        #name = text
        
        if "HOTSPOT" in text:   
            label.set_visible(True)
            label.set_color("red")
        elif "Oncogenic" in text:
            label.set_visible(True)
            label.set_color("blue")
        else:
            label.set_visible(False)
            tick.tick1line.set_visible(False)
            tick.label1.set_visible(False)
            
            
        labels.append(name)
        
    
    r.ax_heatmap.set_xticklabels(labels)
    
    
    
    
    plt.savefig(absoluteFile.replace(".txt", ".pdf"))
    
    
    
    

if __name__ == '__main__':
    main()

