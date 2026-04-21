"""
polyAAPauseScore_FC.py

Modified from polyAAPauseScore.py

Features:
- Input jam format: replicate\tcondition\tjamFile
- Supports gene subset OR ALL genes
- Computes log2FC (Mut/WT)
- Plots per replicate

How to run:
python3 polyAAPauseScore_FC.py seqs.fa jamFiles.txt geneList.txt(or "ALL", for all genes) AA(or "ALL" for all AA)  window size output_prefix

"""

import sys, common, collections, numpy, math
from logJosh import Tee
from pyx import *

print('This is a RAM-hungry script. It will use tens of Gb for C.elegans.')

# -------------------------
# Gene list
# -------------------------

def readGeneList(fn):
    with open(fn) as f:
        return set(line.strip() for line in f if line.strip())

# -------------------------
# Translation
# -------------------------

def convertToProteins(inDict, geneSet):
    aa = {}
    for k, v in inDict.items():
        k2 = k.strip().split()[0]
        if geneSet is None or k2 in geneSet:
            v2 = common.translate(v)[1]
            aa[k2] = v2
    return aa

# -------------------------
# AA window scanning
# -------------------------

def getRuns(protSeqs, N):
    aa = {}
    AAs = ['A','C','D','E','F','G','H','I','K','L',
           'M','N','P','Q','R','S','T','V','W','Y']
    for k, v in protSeqs.items():
        aa[k] = {}
        for ii in range(len(v) - N):
            aa[k][ii] = {}
            temp = v[ii:ii+N]
            for AA in AAs:
                aa[k][ii][AA] = temp.count(AA)
    return aa

# -------------------------
# JAM parsing
# -------------------------

def getPosition(txtList):
    starts = list(set([int(entry.split(':')[1]) for entry in txtList]))
    stops  = list(set([int(entry.split(':')[2]) for entry in txtList]))
    if len(starts) == 1 and len(stops) == 1:
        return starts[0]
    else:
        return 'na'

def parseReads(jamFile, geneSet):
    aa = collections.defaultdict(lambda: collections.defaultdict(int))
    with open(jamFile, 'r') as f:
        f.readline()
        for line in f:
            line = line.strip().split('\t')
            if line[8].endswith(':S'):
                txtList = line[-1].split('|')
                position = getPosition(txtList)
                if position != 'na':
                    for entry in txtList:
                        txtID = entry.split(':')[0]
                        if geneSet is None or txtID in geneSet:
                            aa[txtID][position] += 1
    return aa

def getReads(jamFiles, geneSet):
    aa = []
    with open(jamFiles, 'r') as f:
        for line in f:
            rep, cond, jf = line.strip().split('\t')
            aa.append((rep, cond, parseReads(jf, geneSet)))
    return aa

# -------------------------
# Pair replicates
# -------------------------

def pairReplicates(readTuples):
    pairs = collections.defaultdict(dict)
    for rep, cond, rd in readTuples:
        pairs[rep][cond] = rd

    out = []
    for rep, d in pairs.items():
        if 'WT' in d and 'MUT' in d:
            out.append((rep, d['WT'], d['MUT']))
        else:
            print(f"Warning: replicate {rep} missing WT or MUT")
    return out

# -------------------------
# Metagene (single lib)
# -------------------------

def computeMetaForOne(runDict, readDict, ii, AA, X, N):
    aa = collections.defaultdict(float)
    count = 0

    for txtID, ctDictPerPosition in runDict.items():
        if txtID in readDict:
            for position, ctDictPerAA in ctDictPerPosition.items():
                position2 = position * 3
                if ctDictPerAA[AA] >= ii:
                    tot = 0.0
                    for jj in range(position2 - X, position2 + X + N):
                        tot += readDict[txtID][jj]

                    if tot > 0:
                        count += 1
                        for jj in range(position2 - X - 12, position2 + X + N):
                            aa[jj - position2] += readDict[txtID][jj] / tot

    if count == 0:
        return {}

    cc = {k: v / count for k, v in aa.items()}
    theMetaTotal = sum(cc.values())

    if theMetaTotal == 0:
        return {}

    dd = {k + 12: v / theMetaTotal for k, v in cc.items()}
    return dd

# -------------------------
# Compute log2FC
# -------------------------

def getMetaDataLog2FC(ii, AA, runDict, pairedReads, X, N):
    bb = []
    pseudocount = 1e-6

    for rep, wtDict, mutDict in pairedReads:
        wtMeta  = computeMetaForOne(runDict, wtDict, ii, AA, X, N)
        mutMeta = computeMetaForOne(runDict, mutDict, ii, AA, X, N)

        log2fc = {}
        for pos in range(-X, X + N):
            wt  = wtMeta.get(pos, 0)
            mut = mutMeta.get(pos, 0)
            log2fc[pos] = math.log2((mut + pseudocount) / (wt + pseudocount))

        bb.append((rep, log2fc))

    return bb

# -------------------------
# Prep for plotting
# -------------------------

def getMetaDataPrepared(metaData, X, N):
    aa = []
    for name, metaDict in metaData:
        temp = []
        for ii in range(-X, X + N):
            temp.append((ii, metaDict.get(ii, 0)))
        aa.append((name, temp))
    return aa

# -------------------------
# Plotting
# -------------------------

def mkPlot(ii, AA, runDict, pairedReads, N):
    X = 100

    metaData = getMetaDataLog2FC(ii, AA, runDict, pairedReads, X, N)
    metaDataPrepped = getMetaDataPrepared(metaData, X, N)

    axis_title = f"Pos Rel Start of >= {ii} of {N} are {AA}"

    g = graph.graphxy(
        width=16, height=2, ypos=ii * 4,
        key=graph.key.key(pos='tr', hinside=0),
        x=graph.axis.linear(min=-X, max=X + N, title=axis_title),
        y=graph.axis.linear(min=-2, max=2, title='log2FC (Mut/WT)')
    )

    # zero line
    zero_line = [(-X, 0), (X + N, 0)]
    g.plot(
        graph.data.points(zero_line, x=1, y=2),
        [graph.style.line([style.linestyle.dashed, style.linewidth.Thin])]
    )

    ct = 0
    for name, theData in metaDataPrepped:
        ct += 1
        g.plot(
            graph.data.points(theData, x=1, y=2, title=name),
            [graph.style.line([style.linestyle.solid, common.colors(ct)])]
        )

    return g

# -------------------------
# Main
# -------------------------

def main(args):
    seqFile, jamFiles, geneFile, AAinput, N, outPrefix = args
    N = int(N)

    if geneFile == "ALL":
        geneSet = None
    else:
        geneSet = readGeneList(geneFile)

    cdsSeqs = common.parseFasta(seqFile)
    protSeqs = convertToProteins(cdsSeqs, geneSet)

    runDict = getRuns(protSeqs, N)

    readTuples = getReads(jamFiles, geneSet)
    pairedReads = pairReplicates(readTuples)

    allAAs = ['A','C','D','E','F','G','H','I','K','L',
          'M','N','P','Q','R','S','T','V','W','Y']

    if AAinput == "ALL":
        AAs = allAAs
    else:
        if AAinput not in allAAs:
            print(f"Error: {AAinput} is not a valid amino acid")
            sys.exit()
        AAs = [AAinput]

    for AA in AAs:
        print(f'Working on {AA}')
        c = canvas.canvas()
        for ii in range(3, 7):
            print(f'  Run length {ii}')
            g = mkPlot(ii, AA, runDict, pairedReads, N)
            c.insert(g)

        c.writePDFfile(f'{outPrefix}.{AA}.log2FC')

# -------------------------

if __name__ == '__main__':
    Tee()
    main(sys.argv[1:])
