'''
Author: Juber Patel

'''


import os 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import date
from matplotlib import rcParams
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
from matplotlib.ticker import MaxNLocator
import matplotlib.ticker as mticker






def main():

    d = "/Users/patelj1/current/PromiseStudyJillian/analysis/baseline-8-weeks-ctdna-fraction"

    os.chdir(d)

    tissueFile = "tissue.tsv"
    cfDNAFile = "cfDNA.tsv"

    cfDNA = pd.read_csv(cfDNAFile, sep="\t", index_col=False, keep_default_na=False, low_memory=False)

    cfDNA["PM"] = cfDNA["PatientMutation"].apply(lambda x: x.split(":")[0] + ":" + x.split(":")[1])

    #cfDNA.loc[cfDNA["AF"]>0.05, "AF"] = 0.05

    cfDNA["Waltz_total_t_alt_count"] = cfDNA["Waltz_total_t_alt_count"].astype(int)
    cfDNA.loc[cfDNA["Waltz_total_t_alt_count"]==1, "Waltz_total_t_alt_count"] = 0
    cfDNA.loc[cfDNA["Waltz_total_t_alt_count"]==1, "AF"] = 0

    #cfDNA.loc[cfDNA["Waltz_total_t_alt_count"]>10, "Waltz_total_t_alt_count"] = 10

    
    rcParams['figure.figsize'] = 40,10
    sns.set_theme(style="whitegrid")
    plt.figure(dpi=300)
    plt.tight_layout()

    g = sns.scatterplot(data=cfDNA, x="PM", y="AF", hue="Timepoint")
    plt.xticks(rotation=90)

    plt.tight_layout()
    plt.savefig("IMPACT-mutations-cfDNA-AF.pdf")
    plt.close()

    plt.gca().xaxis.set_major_locator(mticker.MultipleLocator(1))
    g = sns.scatterplot(data=cfDNA, x="PM", y="Waltz_total_t_alt_count", hue="Timepoint")
    plt.xticks(rotation=90)

    plt.gca().xaxis.set_major_locator(mticker.MultipleLocator(1))
    #ax = plt.figure().gca()
    #ax.yaxis.set_major_locator(MaxNLocator(integer=True))

    plt.tight_layout()
    plt.gca().xaxis.set_major_locator(mticker.MultipleLocator(1))
    plt.savefig("IMPACT-mutations-cfDNA-AltCount.pdf")
    plt.close()









if __name__ == '__main__':
    main()
