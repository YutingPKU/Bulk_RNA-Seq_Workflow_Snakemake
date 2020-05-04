#!/usr/bin/python2.7
"""Uppercases the mgi symbol in the csv annotation file"""
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
    #skip the header
    print(f.readline().strip())
    for l in f:
        tmp = l.strip().split(",")
        #uppercase the first element
        tmp[0] = tmp[0].upper()
        print ",".join(tmp)
        #sys.exit()

if __name__=='__main__':
    main()

