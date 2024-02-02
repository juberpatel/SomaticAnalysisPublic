import os 
import sys
import pysam




def main():

    d = "/Users/patelj1/current/RNAMediatedDNARepair/IMPACT"

    os.chdir(d)

    widFile = "t1.tsv"
    MHFile = "mh-wid.tsv"
    
    fastaFile = "/Users/patelj1/resources/impact-GRCh37/Homo_sapiens_assembly19.fasta"

    fasta = pysam.FastaFile(fastaFile)

    #seq = fasta.fetch(region="1:16265786-16265793")
    #print(seq)

    # read wid file line by line

    f = open(widFile)
    o = open(MHFile, "w")
    count = 0
    mh_deletions = 0

    for line in f:
        count += 1

        line = line.rstrip()

        if count == 1:
            line = "\t".join([line, "Left-outer", "Left-inner", "Right-inner", "Right-outer"])
            o.write(line + "\n")
            continue

        words = line.split("\t")
        chr = words[4]
        start = int(words[5])
        end = int(words[6])

        
        reg1 = chr + ":" + str(start) + "-" + str(start+2)
        reg2 = chr + ":" + str(start-3) + "-" + str(start-1)
        reg3 = chr + ":" + str(end-2) + "-" + str(end)
        reg4 = chr + ":" + str(end+1) + "-" + str(end+3)
        

        '''
        reg1 = chr + ":" + str(start) + "-" + str(start+1)
        reg2 = chr + ":" + str(start-2) + "-" + str(start-1)
        reg3 = chr + ":" + str(end-1) + "-" + str(end)
        reg4 = chr + ":" + str(end+1) + "-" + str(end+2)
        '''

        '''
        reg1 = chr + ":" + str(start) + "-" + str(start)
        reg2 = chr + ":" + str(start-1) + "-" + str(start-1)
        reg3 = chr + ":" + str(end) + "-" + str(end)
        reg4 = chr + ":" + str(end+1) + "-" + str(end+1)
        '''

        #print(reg1)
        left_inner = fasta.fetch(region=reg1)
        left_outer = fasta.fetch(region=reg2)
        right_inner = fasta.fetch(region=reg3)
        right_outer = fasta.fetch(region=reg4)

        if left_inner == right_outer or left_outer == right_inner:
            
            line = "\t".join([line, left_outer, left_inner, right_inner, right_outer])
            print(line)
            o.write(line + "\n")
            mh_deletions += 1


    f.close()
    o.close()

    print(f"MH deletions: {mh_deletions}")

        













if __name__ == '__main__':
    main()


