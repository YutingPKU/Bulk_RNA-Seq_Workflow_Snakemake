#!/usr/bin/python2.7
"""converts the ensemble gtf file where the chromosome is named 
10 instead of chr10, or X instead of chrX!"""
import os
import sys
from optparse import OptionParser

def main():
    usage = "USAGE: %prog [options] file1 file2 ... fileN"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-f", "--file", help="file to convert")

    (options, args) = optparser.parse_args(sys.argv)
    
    if not options.file:
        optparser.print_help()
        sys.exit(-1)
    
    f = open(options.file)
    for l in f:
        if l.startswith("#"):
            print(l.strip())
        else:
            tmp = l.strip().split("\t")
            #convert to new chr
            tmp[0] = "chr%s" % tmp[0]
            print "\t".join(tmp)

if __name__=='__main__':
    main()

