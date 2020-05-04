#!/usr/bin/env python
"""Script to convert STAR junction files (SJ.out.tab) to TopHat junctions.bed
format **SO it can be read into IGV.

NOTE: SJ.out.tab format
ref: https://www.biostars.org/p/146700/
column 1: chromosome
column 2: first base of the intron (1-based)
column 3: last base of the intron (1-based)
column 4: strand (0: undefined, 1: +, 2: -)
column 5: intron motif: 0: non-canonical; 1: GT/AG, 2: CT/AC, 3: GC/AG, 4: CT/GC, 5:AT/AC, 6: GT/AT
column 6: 0: unannotated, 1: annotated (only if splice junctions database is used)
column 7: number of uniquely mapping reads crossing the junction
column 8: number of multi-mapping reads crossing the junction
column 9: maximum spliced alignment overhang


NOTE: TopHat junctions.bed format
[seqname] [start] [end] [id] [score] [strand] [thickStart] [thickEnd] [r,g,b] [block_count] [block_sizes] [block_locations]
"start" is the start position of the leftmost read that contains the junction.
"end" is the end position of the rightmost read that contains the junction.
"id" is the junctions id, e.g. JUNC0001
"score" is the number of reads that contain the junction.
"strand" is either + or -.
"thickStart" and "thickEnd" don't seem to have any effect on display for a junctions track. TopHat sets them as equal to start and end respectively.
"r","g" and "b" are the red, green, and blue values. They affect the colour of the display.
"block_count", "block_sizes" and "block_locations":
The block_count will always be 2. The two blocks specify the regions on either side of the junction. "block_sizes" tells you how large each region is, and "block_locations" tells you, relative to the "start" being 0, where the two blocks occur. Therefore, the first block_location will always be zero.

[read_start][junction][read_end]
[block1 ][ ][block2]

"""
import os
import sys
from optparse import OptionParser

def getStrand(s):
    """STAR reports 0 - undefined, 1 = +, 2 = -"""
    if s == "1":
        return "+"
    elif s == "2":
        return "-"
    else:
        return "."

def getLocation(start, end, size):
    """Calculation from https://gist.github.com/fabiolib/ffb21853e3eb3780150074aeffc54901 is:
    END - START + SIZE + 1
    """
    return int(end) - int(start) + int(size) + 1

def main():
    usage = "USAGE: %prog [options] -f [SJ.out.tab file]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-f", "--file", help="STAR SJ.out.tab file to convert")

    (options, args) = optparser.parse_args(sys.argv)
    if not options.file:
        optparser.print_help()
        sys.exit(-1)

    f = open(options.file)
    for (i,l) in enumerate(f):
        tmp = l.strip().split("\t")
        
        #TRANSLATE/MAP
        seqname = tmp[0]
        #start = tmp[1]
        #end = tmp[2]
        start = "%s" % (int(tmp[1]) - int(tmp[8]) - 1)
        end = "%s" % (int(tmp[2]) + int(tmp[8]))
        iid = "JUNC{:07d}".format(i+1)
        score = tmp[7]
        strand = getStrand(tmp[3])
        thickStart = start
        thickEnd = end
        rgb = "255,0,0"
        #BLOCK count ALWAYS 2
        block_count = "2"
        #For blocksizes, I googled "convert sj.out.tab to bed" and found
        #a script that did it this way
        #ref: https://gist.github.com/fabiolib/ffb21853e3eb3780150074aeffc54901
        block_sizes="%s,%s" % (tmp[8],tmp[8])
        block_loc="0,%s" % getLocation(tmp[1], tmp[2], tmp[8])

        print("\t".join([seqname, start, end, iid, score, strand, thickStart,
                         thickEnd, rgb, block_count, block_sizes, block_loc]))

if __name__=='__main__':
    main()

