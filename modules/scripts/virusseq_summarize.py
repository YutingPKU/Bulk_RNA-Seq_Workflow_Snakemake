#!/usr/bin/env python
"""Script to summarize the virusseq results--outputs to stdout: AS CSV
SampleID, TranscriptID, Counts, FPKM
"""

import os
import sys
from optparse import OptionParser

def main():
    usage = "USAGE: %prog -f "
    optparser = OptionParser(usage=usage)
    optparser.add_option("-f", "--file", help="virusseq.table file")
    (options, args) = optparser.parse_args(sys.argv)

    if not options.file:
        optparser.print_help()
        sys.exit(-1)

    f=open(options.file)
    #BURN
    tmp = f.readline()
    #RELYING on the table file already being sorted!
    lastID = ''
    #print header
    print(",".join(["SampleID","TranscriptID","Counts","FPKM"]))
    for l in f:
        tmp = l.strip().split(',')
        #FPKM - only one sig. digit
        tmp[-1] = "%.1f" % float(tmp[-1])
        if tmp[0] != lastID:
            #NEW SET
            ct = 1
            lastID = tmp[0]
            #PRINT the whole row
            print(",".join(tmp))
            continue

        if ct < 5:
            #PRINT JUST the new info
            print(",".join(tmp))
            tmp[0] = ''
            ct += 1

if __name__=='__main__':
    main()


