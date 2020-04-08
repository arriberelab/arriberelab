import os
OUTPUT_FILE = ".summary.txt"
# import theREALcountReads as jkim
# jkim.main(outPrefix, fastqFile)

def countTrimmed(filename):
    with open(filename, 'r') as a:  # read filename as a
        count = 0
        for line in a:  # increase readCount by 1 for every line in a that starts with '@'
            if line.startswith('@'):
                count += 1
    return count

def countSAMS(filename):
    with open(filename, 'r') as b:
        samcount=0
        for line in b:
            samcount+=1
    return samcount


def main(outPrefix, fastq):
    # try:
    #     os.remove(outPrefix+OUTPUT_FILE) #remove file to rewrite
    # except Exception as e: #if no file exists then exception occurs instead
    #     print(f'ERROR: {e}')
    #     print('Running now')
    input_files = {}
    l = ['.trimmed.collapsed.fastq', '.trimmed.tooShort.fastq', '.trimmed.tooLong.fastq', '.joshSAM'] #list of filenames
    input_files[outPrefix] = {} #dictionary of outPrefixes
    fastqCount = countTrimmed(fastq) #passing fastq through countTrimmed and assign to variable fastq
    input_files[outPrefix]['.fastq'] = fastqCount #save fastqcount to '.fastq' key in outPrefix dictionary
    for i in l:
        input_files[outPrefix][i] = 0 #initialize to 0
        filename = outPrefix+i

        if i == '.joshSAM':
            readCount = countSAMS(filename)
        elif i == '.trimmed.tooLong.fastq' or '.trimmed.tooShort.fastq' or '.trimmed.collapsed.fastq':
            readCount = countTrimmed(filename)
        # else:
        #     print('u thought')
        #     exit()
        input_files[outPrefix][i]= readCount #update with readCount
    #print(input_files)
    #print(input_files[outPrefix])
    with open(outPrefix+OUTPUT_FILE, 'w+') as f:
        for key in input_files.keys():
            f.write(f'{key}: \n')
            for keyprime in input_files[key].keys():
                percentage= round((input_files[key][keyprime] / input_files[key]['.fastq'])* 100, 2)
                f.write(f'{keyprime}:{input_files[key][keyprime]} "reads", {percentage}%\n')


if __name__ == "__main__":
    outPrefix = 'JK0043_15-18mers'
    fastq='JK0043_S38_L006_R1_001.fastq'
    main(outPrefix, fastq)
