__author__ = 'RobbyBoney'
"""
DATE_Started: FEB 2016
DATE_Updated: JAN 9 2017
VERSION: alpha2

OVERVIEW: this script is focussed soley on graphing and visualizng data generated from bio_inforegen2.py

"""
import os
import time
import numpy as np
import matplotlib.pyplot as plt

def ARO_dotplot(domfile, hitfile):
    """
    OVERVIEW:
        this funtion takes 2 files and produces the graph of divergence from hits to domain hits

    USE of FUNC:
        domfile ===> 

    """
    datapoints = []
    data_x = []
    data_y = []

    domfile = open(domfile, "r")
    domfile = domfile.read()
    domfile_s = domfile.splitlines()

    hitfile = open(hitfile, "r")
    hitfile = hitfile.read()
    hitfile_s = hitfile.splitlines()

    for line in hitfile_s:
        
        if (line):
            lnID = line.split()[3][:-4]
            
            for dom in domfile_s:
                dmID = dom.split()
                if (lnID == dmID[5]):
                    #print(lnID, "= ", line.split()[0],"---", dmID[5], "= ", dmID[0])
                    data_x.append(float(line.split()[0]))
                    data_y.append(float(dmID[0]))
                    break

    datapoints.append(data_x)
    datapoints.append(data_y)
    return datapoints


if __name__ == "__main__":
    conhits = "/Users/themusicman/Projects/Python/bioinformaticsJAN2017/Results/Jan17-2017_condensedHits.txt"
    domhits = "/Users/themusicman/Projects/Python/bioinformaticsJAN2017/Results/Feb5-2017_domainHits_con_sorted.txt"

    resultsL = ARO_dotplot(domhits, conhits)
    #resultsL = [ [0.5,0.3,0.4], [0.5,0.3,0.4] ]
    bounds = [0, 1]

    plt.plot(bounds, "k-") 
    plt.axis([0, 0.4, 0, 0.4], )
    plt.scatter(resultsL[0], resultsL[1], s=150, c='b', marker='o', linewidths=0.2)
    
    #plt.tight_layout(pad=3, h_pad=None, w_pad=None, rect=None)
    plt.show()

