'''
Created on Nov 15, 2020

@author: Juber Patel
'''

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import matplotlib.ticker as mtick
import pandas as pd
import os
from numpy import float64
import re


def plotAverages(name, pairs, total, colors):
    
    # y-axis in bold
    #rc('font', weight='bold')
    plt.rcParams.update({'font.size': 15})
    plt.rcParams.update({'font.weight': "bold"})
    plt.rcParams.update({'figure.dpi': 300})
    
    fig, ax = plt.subplots()
    
    bot = 0.0;
    index=0
    for exposure, sig in pairs:
        exposure = float64(exposure)
        exposurePercent = (exposure/total) * 100
        container = ax.bar([""], exposurePercent, width=1, label=sig, bottom=bot, color=colors[sig], edgecolor='white')
        h = container[index].get_height()
        w = container[index].get_width()
        plt.text(container[index].get_x() + w / 2., bot + h / 2., 
                 "%s" % (sig + "(" + "{:.0f}".format(exposure) + ")"), ha="center", 
                 va="center", color="black", fontsize=25, fontweight="bold")
        bot += exposurePercent
    
    ax.set_ylabel('Exposure')
    ax.yaxis.set_major_formatter(mtick.PercentFormatter())
    
    #ax.legend()
    #plt.show()
    
    plt.tight_layout()
    plt.savefig(name + ".pdf")
    plt.close()
    
    
    
    

def main():
    os.chdir("/Users/patelj1/current/ACCESSGRAILMetasaticBreastCancer/analysis/signatures-cfDNA")
    MutationalPatternsFile = "MutationalPatterns.tsv"
    
    
    signatureUniverseFile = "/Users/patelj1/current/workspace/SomaticAnalysis/signature-universe-aggregated.txt"
    
    
    
    # load signature universe and assign colors
    cmap = plt.get_cmap("tab20")
    colors = {}
    counter = 0.0;
    offset = 1.0/40.0
    with open(signatureUniverseFile) as f:
            for line in f:
                words = re.split('\t|\n', line)
                colors[words[0]] = cmap(counter/20.0 + offset)
                counter += 1
                
    # create aggregate signature dictionary
    if signatureUniverseFile.endswith("signature-universe-aggregated.txt"):        
        aggSignatures = {}
        with open(signatureUniverseFile) as f:
            for line in f:
                words = re.split('\t|\n', line)
                sig = words[0]
                if sig == "Other":
                    continue
                contributors = words[1].split(" ")
                for contributor in contributors:
                    aggSignatures["SBS" + str(contributor)] = sig
                    
        print(aggSignatures)
                    
    
    with open(MutationalPatternsFile) as f:
            l = f.readline().strip()
            names = re.split('\t|\n', l)
            
            for line in f:
                s = []
                e = []
                words = re.split('\t|\n', line.strip())
                sample = words[0]
                total = 0.0
                for i in range(1, len(words)):
                    total += float64(words[i])
                    
                if(total < 10):
                    continue
                
                for i in range(1, len(words)):
                    v = float64(words[i])
                    #if(v/total < 0.10):
                        #continue
                    e.append(v)
                    s.append(names[i])
                
                pairs = tuple(zip(float64(e), s))
                pairs = sorted(pairs, reverse=True)
                print(sample)
                print(pairs)
                
                # aggregate signatures
                if signatureUniverseFile.endswith("signature-universe-aggregated.txt"):
                    aggregated = {}
                    total = 0.0
                    for exposure, sig in pairs:
                        total += exposure
                        aggS = aggSignatures.get(sig, "Other")
                        aggE = aggregated.get(aggS, 0.0)
                        aggregated[aggS] = aggE + exposure
                    
                    # anything below threshold is given to 'Other'
                    t = {'Other': 0.0}
                    for aggS, aggE in aggregated.items():
                        if aggE/total < 0.10 or aggS == "Other":
                            t['Other'] = t['Other'] + aggE
                        else:
                            t[aggS] = aggE
                    
                    if t['Other']/total < 0.10:
                        t.pop('Other')
                        
                    aggregated = t
                        
                    # convert dictionary to tuple
                    s = list(aggregated.keys())
                    e = []
                    for aggS in s:
                        e.append(aggregated[aggS])
                    
                    pairs = tuple(zip(e, s))
            
                pairs = sorted(pairs, reverse=True)
                print(pairs)
                plotAverages(sample, pairs, total, colors)




if __name__ == '__main__':
    main()
    
    
    
