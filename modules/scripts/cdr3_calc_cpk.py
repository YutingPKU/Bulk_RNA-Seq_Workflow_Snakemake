#!/usr/bin/env python
"""Script to collect the trust results across all samples. outputs to stdout
SampleID, AssemblyCount, TotalCount, CPK
Where:
AssemblyCount = # of lines in trust .fa (account for headers if any) / 2
**Don't forget to divide by 2!
Total count = 1st number in est_lib_size (4th field of first line)

CPK = (AssemblyCount*1000)/TotalCount
We graph CPK as a box plot and put into viper report

OUTPUTS to stdout
"""

import os
import sys
from optparse import OptionParser

def main():
    usage = "USAGE: %prog -f [trust .fa output FILE_1] -f [trust .fa output FILE_2] ...-f [trust .fa output FILE_N]\n -d [directory of output file DIR_1] -d [directory of output file DIR_2]...."
    optparser = OptionParser(usage=usage)
    optparser.add_option("-f", "--files", action="append", help="list of trust .fa files")
    optparser.add_option("-d", "--dirs", action="append", help="list of trust output directories")
    (options, args) = optparser.parse_args(sys.argv)

    if not options.files and not options.dirs:
        optparser.print_help()
        sys.exit(-1)

    #COLLATE the files and dir files into one big list
    _trustFiles = []
    if options.files:
        _trustFiles = options.files
    if options.dirs:
        for d in options.dirs:
            ls = [os.path.join(d, f) for f in os.listdir(d)]
            _trustFiles.extend(ls)
    
    print(",".join(["SampleID","AssemblyCount","TotalCount","CPK"]))

    #CHECK for valid files
    _trustFiles = list(filter(lambda f: f.endswith(".fa"), _trustFiles))
    for (i, trust_f) in enumerate(_trustFiles):
        #NOTE: handle header
        num_lines = 0
        skip_ct = 0
        totalCount = 0
        sampleID = None
        
        f = open(trust_f)
        for l in f:
            num_lines += 1
            if l.startswith("#"): #skip comment
                skip_ct += 1
            else:
                tmp = l.strip().split("+")                
                #is totalCount set?
                if not totalCount:
                    totalCount = int(tmp[3].split("=")[1].split("-")[0])
                #did we get sampleId
                if not sampleID:
                    firstElm = tmp[0][1:] #skip the > in front
                    sampleID = firstElm.split('.')[0]
        #finish reading the file
        f.close()
        #good reading? --otherwise skip
        if totalCount: 
            assemblyCount = (num_lines - skip_ct) / 2.0
            CPK = (assemblyCount*1000)/totalCount
            print(",".join([sampleID, "%.1f" %assemblyCount, str(totalCount), "%.2f" % CPK]))

if __name__=='__main__':
    main()


