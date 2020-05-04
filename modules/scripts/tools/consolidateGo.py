#!/usr/bin/python2.7
"""GIVEN a biomart output, consolidate the rows to JUST one row per gene
NOTE: the problem is in the GO IDs and Go Terms

TRY: to consolidate them both with a ;

NOTE: gene ids map to multiple ensemblid (773) and entrezid (190)
FOR THESE duplicates, take the first ENSEMBLID
"""
import os
import sys
from optparse import OptionParser

def main():
    usage = "USAGE: %prog [options] -f [FILE]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-f", "--file", help="file to consolidate")

    (options, args) = optparser.parse_args(sys.argv)
    
    if not options.file:
        optparser.print_help()
        sys.exit(-1)
    
    d = {}
    f = open(options.file)
    #skip the header
    print(f.readline().strip())
    for l in f:
        tmp = l.strip().split(",")
        if not tmp[0] in d:
            #STORE the first 3 cols, and the GO IDs and TERMS as list
            d[tmp[0]] = [tmp[0], tmp[1], tmp[2], tmp[3], [tmp[4]], [tmp[5]]]
        else:
            #EXISTING gene--add to the GO LISTS
            d[tmp[0]][4].append(tmp[4])
            d[tmp[0]][5].append(tmp[5])

    #print out the genes as a csv
    for k in d.keys():
        #REMOVE '' from lists
        if '' in d[k][4]:
            d[k][4].remove('')
        if '' in d[k][5]:
            d[k][5].remove('')

        #print set(d[k][4])
        GO_terms = ";".join(set(d[k][4]))
        GO_ids = ";".join(set(d[k][5]))
        print ",".join([d[k][0], d[k][1], d[k][2], d[k][3], GO_terms, GO_ids])

if __name__=='__main__':
    main()

