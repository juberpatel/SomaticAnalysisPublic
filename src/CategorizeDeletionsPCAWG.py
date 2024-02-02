'''
Created on Apr 22, 2022

@author: patelj1
'''

#import portion as p

import os
import sys
import intervaltree as it



def main():
    
    d = "/Users/patelj1/current/RNAMediatedDNARepair/PCAWG/"    
    os.chdir(d)
    
    #exonsFile = "t1.txt"
    exonsFile = "juber-hg19-gene-list.bed"
    deletionsFile = "PCAWG-ICGC-deletions-gt-eq-10bp.tsv"
    outFile = deletionsFile.replace(".tsv", "-categorized.tsv")
    
    genes = set()
    max = {}
    min = {}
    chr = {}
    
    intervalTrees = {}
    geneTrees = {}
    
    # create intervalTrees
    for i in range(1, 23):
        intervalTrees[str(i)] = it.IntervalTree()
        geneTrees[str(i)] = it.IntervalTree()
        
    intervalTrees["X"] = it.IntervalTree()
    intervalTrees["Y"] = it.IntervalTree()
    
    geneTrees["X"] = it.IntervalTree()
    geneTrees["Y"] = it.IntervalTree()
    
    
    
    counter = 0

    # read exons file
    with open(exonsFile) as f:
        for line in f:
            tokens = line.rstrip("\n").split("\t")
            
            if tokens[4].startswith("snp"):
                continue
            
            if tokens[4].startswith("MIR"):
                continue
            
            if "exon" not in tokens[4]:
                continue
            
            if tokens[0] not in intervalTrees:
                continue
            
            # add the gene to the gene set
            gene = tokens[4].split("_")[0]
            genes.add(gene)
            chr[gene] = tokens[0]
            
            start = int(tokens[1])
            end = int(tokens[2]) + 1
            intervalTrees[tokens[0]][start:end] = tokens[4]
            
            # update gene range
            if gene not in min or start < min[gene]:
                min[gene] = start
                
            if gene not in max or end > max[gene]:
                max[gene] = end
            
            counter += 1
            print(counter)
            
    
    f.close()
    
    
    # build gene trees
    for gene in min:
        c = chr[gene]
        geneTrees[c][min[gene]:max[gene]] = gene
        
    
    
    '''
    #13    26975485    26975602
    c = "13"
    position = 26984149
    overlap = sorted(geneTrees[c][position])
    
    print(overlap)
    sys.exit()
    '''
    
    
    
    
    # open output file
    w = open(outFile, "w")
    
    counter = 0
    # read the deletions file
    with open(deletionsFile) as f:
        for line in f:
            
            if counter == 0:
                w.write(line.rstrip("\n") + "\tDeletion_Category\n")
                counter += 1
                continue
            
            
            category = "UNKNOWN"
            tokens = line.rstrip("\n").split("\t")
            t = intervalTrees[tokens[1]]
            
            # intron deletion
            found = False
            for i in range(-2, 3):
                if found:
                    break
                
                for j in range(-2, 3):
                    start = int(tokens[2]) + i
                    end = int(tokens[3]) - j + 1
                    overlap = sorted(t[start:end])
                    #print(str(start) + "\t" + str(end) + "\t" + str(overlap))
                    
                    if len(overlap) != 2:
                        continue
                    
                    interval = overlap[0]
                    interval1 = overlap[1]
                    if abs(int(tokens[2])-interval.end) <= 2 and abs(int(tokens[3])-interval1.begin) <= 2:
                        category = "intron deletion"
                        w.write(line.rstrip("\n") + "\t" + category + "\n")
                        counter += 1
                        found = True
                        break
            if found:
                continue
                    
                    
            
            # exon deletion
            found = False
            for i in range(-2, 3):
                if found:
                    break
                
                for j in range(-2, 3):
                    start = int(tokens[2]) + i
                    end = int(tokens[3]) - j + 1
                    overlap = sorted(t[start:end])
                    #print(str(start) + "\t" + str(end) + "\t" + str(overlap))
                    
                    if len(overlap) != 1:
                        continue
                    
                    interval = overlap[0]
                    if abs(int(tokens[2])-interval.begin) <= 2 and abs(int(tokens[3])-interval.end) <= 2:
                        category = "exon deletion"
                        w.write(line.rstrip("\n") + "\t" + category + "\n")
                        counter += 1
                        found = True
                        break
            if found:
                continue
            
            
            # genes absent in the annotation: assume exonic variants
            if tokens[0] not in genes:
                category = "exonic"
                w.write(line.rstrip("\n") + "\t" + category + "\n")
                counter += 1
                continue
                
            
            start = int(tokens[2])
            end = int(tokens[3]) + 1
            overlap = sorted(t[start:end])
            
             
            # if the overlap is empty, it's an intronic or intergenic deletion
            if len(overlap) == 0:
                
                # decide if intronic or intergenic
                p = sorted(geneTrees[tokens[1]][start])
                
                if len(p) == 0:
                    category = "intergenic"
                else:
                    category = "intronic"
                
                w.write(line.rstrip("\n") + "\t" + category + "\n")
                counter += 1
                continue
                
            
            interval = overlap[0]
            # if contained within the first interval, it's an exonic deletion
            if interval.begin < start and interval.end > end:
                category = "exonic"
                w.write(line.rstrip("\n") + "\t" + category + "\n")
                counter += 1
                continue
            
            # if bigger than first interval, it is intron-exon spanning
            if start <= interval.begin or end >= interval.end:
                category = "intron-exon spanning"
                w.write(line.rstrip("\n") + "\t" + category + "\n")
                counter += 1
                continue
            
            
               
                
            #print(tokens[0] + "\t" + tokens[4] + "\t" + tokens[5] + "\t" + tokens[6])
            #print(overlap)
            
            
                
            
    f.close()
    w.close()
    
    
    
    
    
    
    
    



if __name__ == '__main__':
    main()