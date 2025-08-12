"""
converted to python 3 on 4/17/2020

Joshua Arribere Oct 7, 2013
Marcus Viscardi Jan 9, 2024 - At some point someone updated this to work with jam files rather than JoshSam

Script to get gene counts for both Sense and Antisense-mapping reads separately, per gene.
    Will output gene ct file, and then run DESeq to median normalize the same file.

Input: libsFile.txt - tab-delimited file of the format
        name1   cond1    file1.jam
        ...
        namen   condj   filen.jam
    Here condi is used to combine all libraries from the same condition. e.g. L1 libraries or mutant libraries.
        condi can be any string. It will be used for DESeq, whereas namei will appear as the header in the outFile.
    relStart,relStop - two numbers that will be used to cutoff any reads that fall outside of these bounds.
        relStart is the 5'position of the read relative to the Start codon, and relStop same for stop codon. Because
        I do gene cts across the entire gene, which may not contain one single CDS, I use the intersection of all gene's
        txts' CDSs for the bounds.

Output: geneCts - files for each of sense/antisense
    sizeFactors - From DESeq, estimates of corrective factors for library size differences
    geneCtsMedianNormalized - geneCts/sizeFactors for each gene

run as python medianNormalizerWithDESeq.py [options] libsFile.txt outPrefix
EDIT: Oct 8, 2013 - JOSH added to include CDS bound option
EDIT: March 10, 2014 - JOSH edited to preserve library input order
"""
import sys, common, parse_joshSAM, csv, os
from logJosh import Tee
from optparse import OptionParser
import numpy, collections

def writeGeneCtFile(outPrefix,geneCtDict,allgenes):
    """geneCtDict={header:[cond,{gene:ct}]}. Will write a tab-delimited file of gene cts. allgenes is a dict
    of all genes used in geneCtDict"""
    
    with open(outPrefix+'.geneCt','w') as f:
        fwriter=csv.writer(f,delimiter='\t')
        headers=[entry[0] for entry in geneCtDict]#changed March 20, 2014
        fwriter.writerow(['']+headers)
        for gene in allgenes:
            theRow=[gene]
            for ii in range(len(headers)):#use headers, not geneCtDict here b/c headers is guaranteed to be in the same order as when I wrote the top of the file
                try:#this is currently unnecessary b/c geneCtDict[header][1] is a defaultdict
                    theRow.append(geneCtDict[ii][2][gene])
                except KeyError:
                    #print 'ehllo'
                    theRow.append(0.)
            fwriter.writerow(theRow)
    
    #now write the conditions file for DESeq
    with open(outPrefix+'.conditions','w') as f:
        fwriter=csv.writer(f,delimiter='\t')
        fwriter.writerow([geneCtDict[ii][1] for ii in range(len(headers))])
    
    #now write the libTypes file for DESeq
    with open(outPrefix+'.libTypes','w') as f:
        fwriter=csv.writer(f,delimiter='\t')
        fwriter.writerow(['single-end']*len(headers))

def getGeneCts(libsDict,outPrefix,SorAS,bounds=False,diffExp=False):
    """Will get gene counts for SorAS-strand mapping reads for each dict in libsDict={name:[cond,file]}"""
    
    #first get the gene cts.
    aa=[]#changed March 20, 2014
    bb=[]
    genes={}#this will keep track of all genes
    for entry in libsDict:
        name=entry[0]#changed March 20, 2014
        print(name, entry[1])#changed March 20, 2014
        libFile=entry[2]#changed March 20, 2014
        if SorAS=='S':
            geneCts=parse_joshSAM.getGeneReadCts(libFile,senseOnly=True,cdsBounds=bounds)
        elif SorAS=='AS':
            geneCts=parse_joshSAM.getGeneReadCts(libFile,antisenseOnly=True,cdsBounds=bounds)
        total=sum(geneCts.values())/1000000.
        aa.append([name,entry[1],geneCts])#changed March 20, 2014
        bb.append([name,entry[1],dict((gene,geneCts[gene]/total) for gene in geneCts)])
        for gene in geneCts:
            genes[gene]=1
        
        total=sum(geneCts.values())
        print(total)
    
    outPrefix2=outPrefix+'_'+SorAS
    #Now write the first outFile
    print('Writing gene count file to '+outPrefix2)
    #print 'Writing gene count file as RPM, NOT raw counts. to change, switch bb to aa in below line'
    writeGeneCtFile(outPrefix2,aa,genes)
    
    #Now run DESeq
<<<<<<< HEAD
    #rscriptLocation='/data15/joshua/github/200329_arribereLabPipeline/step2_analysisScripts/'
    rscriptLocation= '/data14/chloe/scripts/arribere_github/arriberelab/step2_analysisScripts/'
=======
    # Added by Marcus on 1/9/2024, this will look for
    #   the median_normalize_DESeq2.r script in the same
    #   directory as this script!!
    rscriptLocation = os.path.dirname(os.path.realpath(__file__))
>>>>>>> c082f6a2d514550d84c6d94d8a6ec29594df7fca
    if not diffExp:
        os.system(f'Rscript {rscriptLocation}median_normalize_DESeq2.r '+outPrefix2+'.geneCt '+\
                  outPrefix2+'.conditions '+outPrefix2+'.libTypes '+outPrefix2+'.DESeqgeneCts')
    else:
        os.system(f'Rscript {rscriptLocation}median_normalize_DESeq2.r '+outPrefix2+'.geneCt '+\
                  outPrefix2+'.conditions '+outPrefix2+'.libTypes '+outPrefix2+'.DESeqgeneCts 1')
    
    #Remove the useless intermediate files
    #os.system('rm '+outPrefix2+'.conditions')
    #os.system('rm '+outPrefix2+'.libTypes')

def main(args):
    #parse cdsBounds (as an option)
    parser = OptionParser()
    parser.add_option("-b", "--cdsBounds", type="int", nargs=2, dest="cdsBounds",
                      help="will look only at reads within these bounds")
    
    #added March 4, 2014
    parser.add_option("-d", "--diffExpr", type=int, nargs=1, dest="diffExpr",default=0,
                      help="if set to True, will run differential expression analysis")
    
    (options, args) = parser.parse_args()
    cdsBounds=options.cdsBounds
    diffExpr=options.diffExpr
    print(cdsBounds, ' cdsBounds are inclusive.')
    libsFile,outPrefix=args[0:]
    
    #Rough parsing of the libsFile
    libsDict=[]#changed March 20, 2014
    with open(libsFile,'r') as f:
        for line in f:
            line=line.strip().split()
            libsDict.append([line[0],line[1],line[2]])#changed March 20, 2014
    
    if cdsBounds:
        outPrefix=outPrefix+'_Bounds_'+'_'.join(map(str,cdsBounds))
    print('Writing to: '+outPrefix)
    
    #First do the sense
    getGeneCts(libsDict,outPrefix,'S',bounds=cdsBounds,diffExp=diffExpr)
    
    #Next do the antisense
    getGeneCts(libsDict,outPrefix,'AS',bounds=cdsBounds,diffExp=diffExpr)

if __name__=='__main__':
    Tee()
    main(sys.argv[1:])
