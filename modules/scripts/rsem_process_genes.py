#!/usr/bin/env python
"""Script to process RSEM genes.results file b/c the gene id is duplicated
"""

import os
import sys
from optparse import OptionParser

def main():
    usage = "USAGE: %prog -f rsem.genes.results"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-f", "--file", help="file to process")
    (options, args) = optparser.parse_args(sys.argv)

    if not options.file:
        optparser.print_help()
        sys.exit(-1)

    f = open(options.file)
    #Handle the header
    l = f.readline().strip()
    print(l)
    for l in f:
        tmp = l.strip().split('\t')
        gene_id = tmp[0].split("_")
        #REPLACE the duplicated gene ids w/ singletons
        tmp[0] = gene_id[0]
        print("\t".join(tmp))

if __name__=='__main__':
    main()


