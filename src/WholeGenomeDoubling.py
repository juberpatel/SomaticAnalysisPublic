'''
Created on Jun 22, 2022

@author: Juber Patel
'''

import pandas as pd
import os
import seaborn as sb
import copy
import matplotlib.pyplot as plt

def main():
    
    d = "/Users/patelj1/current/Shirin-Metaplastic"
    
    os.chdir(d)
    
    f = "metaplastic-floh-ploidy.txt"
    
    df = pd.read_csv(f, sep="\t", index_col=False, keep_default_na=False, low_memory=False)
    
    # add WGD verdict
    # Equation of the line: 2X + Y = 2.9
    df["WGDValue"] = df["fLOH"]*2 + df["Ploidy"]
    df["WGD"] = "No WGD"
    df.loc[df["WGDValue"] > 2.9, "WGD"] = "WGD"
    
    df.to_csv(f.replace(".txt", "-wgd.txt"), sep="\t", index=False)
    
    
    
    fig = plt.figure(dpi=300)
    plt.tight_layout()
    
    
    sb.set_style("whitegrid")
    
    p = sb.scatterplot(data=df, x="fLOH", y="Ploidy", hue="WGD", alpha=0.7)
    p.set_xlabel("Fraction of genome with LOH")
    p.set(title='Metaplastic (21 samples)')
    
    plt.plot([0, 1.45], [3, 0],  "--", linewidth=2, color="black")
    plt.xlim(0, 1)
    plt.ylim(0, 10)
    
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)
    plt.tight_layout()
    
    plt.savefig(f.replace(".txt", ".pdf"))
    
    
    



if __name__ == '__main__':
    main()
