#!/usr/bin/env python
"""Script to convert mouse gene names to their human ortholog
*REQUIRES a mapping csv file, where the first col are mouse genes and the 2nd
col is the human ortholog

Input: takes a csv file, and assumes the 0th col is the mouse gene col
OUTPUT: will replace the mouse gene name with the human equivalent IF there is
        a human ortholog, otherwise it will leave it in place
"""

import os
import sys
from optparse import OptionParser

def main():
    usage = "USAGE: %prog -m [map file] -i [input csv file- first col are mouse genes]  -o [output file name]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-m", "--map", help="map file")
    optparser.add_option("-i", "--input", help="input file")
    optparser.add_option("-o", "--out", help="output file")
    (options, args) = optparser.parse_args(sys.argv)

    if not options.map or not options.input or not options.out:
        optparser.print_help()
        sys.exit(-1)

    #READ in the map
    #NOTE: uppercaseing the gene names so there's no mismatch there
    f = open(options.map)
    hdr = f.readline()
    nameMap = {}
    for l in f:
        tmp = l.strip().split(",")
        nameMap[tmp[0].upper()] = tmp[1].upper()
    f.close()
    
    out = open(options.out, "w")
    #READ and replace the names-- ensuring there are no duplicates
    #b/c multiple mm genes map to a single hs gene
    in_f = open(options.input)
    duplicate_ls = []
    for l in in_f:
        tmp = l.strip().split(",")
        mm = tmp[0].upper()
        if mm in nameMap:
            tmp[0] = nameMap[mm].upper()
        #check for duplicates
        if tmp[0] not in duplicate_ls:
            out.write("%s\n" % ",".join(tmp))
            duplicate_ls.append(tmp[0])
    in_f.close()
    out.close()

if __name__=='__main__':
    main()


