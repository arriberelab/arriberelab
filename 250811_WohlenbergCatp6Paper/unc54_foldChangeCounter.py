"""
Modified by CW on 240611 to be customizable to any wt and mut lis
Script to calculate unc-54 fold change relative to wildtype lib.
Normalizes by six RPGs.

Run as python3 unc54_foldChangeCounter.py wtLib.jam mutantLib.jam outPrefix
"""

import os,sys,csv,statistics,subprocess

def main(args):
    wtLib,mutLib,outPrefix=args[0:]
    
    print(f"wtLib {wtLib}")
    wtUnc = os.system(f"grep -c 'unc-54' {wtLib}")
    wtRps6 = os.system(f"grep -c 'WBGene00004475' {wtLib}")
    wtRps12 = os.system(f"grep -c 'WBGene00004481' {wtLib}")
    wtRps24 = os.system(f"grep -c 'WBGene00004493' {wtLib}")
    wtRpl6 = os.system(f"grep -c 'WBGene00004417' {wtLib}")
    wtRpl12 = os.system(f"grep -c 'WBGene00004424' {wtLib}")
    wtRpl24 = os.system(f"grep -c 'WBGene00004436' {wtLib}")
    
    print(f"mutLib {mutLib}")
    mutUnc = os.system(f"grep -c 'unc-54' {mutLib}")
    mutRps6 = os.system(f"grep -c 'WBGene00004475' {mutLib}")
    mutRps12 = os.system(f"grep -c 'WBGene00004481' {mutLib}")
    mutRps24 = os.system(f"grep -c 'WBGene00004493' {mutLib}")
    mutRpl6 = os.system(f"grep -c 'WBGene00004417' {mutLib}")
    mutRpl12 = os.system(f"grep -c 'WBGene00004424' {mutLib}")
    mutRpl24 = os.system(f"grep -c 'WBGene00004436' {mutLib}")
    
    allRPGs=[]
    rps6Norm=wtRps6/mutRps6
    allRPGs.append(rps6Norm)
    rps12Norm=wtRps12/mutRps12
    allRPGs.append(rps12Norm)
    rps24Norm=wtRps24/mutRps24
    allRPGs.append(rps24Norm)
    rpl6Norm=wtRpl6/mutRpl6
    allRPGs.append(rpl6Norm)
    rpl12Norm=wtRpl12/mutRpl12
    allRPGs.append(rpl12Norm)
    rpl24Norm=wtRpl24/mutRpl24
    allRPGs.append(rps24Norm)
    
    avgNorm=statistics.mean(allRPGs)
    mutUncNorm=mutUnc*avgNorm
    mutUncFoldChg=mutUncNorm/wtUnc
    with open (outPrefix + '.txt','w') as out:
        outwriter=csv.writer(out,delimiter='\t')
        outwriter.writerow([wtUnc,mutUnc,avgNorm,mutUncNorm,mutUncFoldChg])
if __name__ == '__main__':
    main(sys.argv[1:])
