#!/usr/bin/env python
"""Script to collect the FPKM results (isoforms.fpkm_tracking, and maybe 
later genes.fpkm_tracking) across all samples
Outputs to stdout:
(Transcript/Gene)ID, Sample1, ..., SampleN
XXX_TRANSCRIPT_ID, SAMPLE1_FPKM, ..., SampleN_FPKM
"""

import os
import sys
from optparse import OptionParser

def main():
    usage = "USAGE: %prog -f [FPKM FILE1] -f [FPKM FILE2] ... -f [FPKM FILE N]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-f", "--fpkms", action="append", help="list of isoform.fpkm_tracking or gene.fpkm_tracking cufflink files")
    optparser.add_option("-n", "--name", default="ID", help="name to call the first col, e.g. Gene_ID or Transcript_ID (default: ID)")
    optparser.add_option("-i", "--id_col", default=0, help="Column in the files that correspond to the id (default: 0)")
    optparser.add_option("-c", "--fpkm_col", default=9, help="Column in the files that correspond to the FPKM (default: 9)")

    (options, args) = optparser.parse_args(sys.argv)

    if not options.fpkms:
        optparser.print_help()
        sys.exit(-1)

    #rename:
    iid = int(options.id_col)
    fid = int(options.fpkm_col)

    matrix = {}
    fpkm_f = options.fpkms[0]
    #INIT the samples list w/ the first sample name
    samples = [fpkm_f.split("/")[-1].split(".")[0]]
    fpkm_f = open(fpkm_f)
    #BURN: read the header
    l = fpkm_f.readline()
    for l in fpkm_f:
        tmp = l.strip().split("\t")
        #initialize: matrix to {ID: [FPKM_1, 0, 0, ..., 0]}
        matrix[tmp[iid]] = [tmp[fid]] + [0]*(len(options.fpkms)-1)
    fpkm_f.close()

    #DEAL with rest:
    for (i, fpkm_f) in enumerate(options.fpkms[1:]):
        #APPEND to the samples name
        samples.append(fpkm_f.split("/")[-1].split(".")[0])
        fpkm_f = open(fpkm_f)
        #BURN: read the header
        l = fpkm_f.readline()
        for l in fpkm_f:
            tmp = l.strip().split("\t")
            if (tmp[iid] in matrix):
                #ADD value
                matrix[tmp[iid]][i+1] = tmp[fid]
            else:
                #NEW element--does this ever happen??
                foo = [0]*len(options.fpkms)
                foo[i+1] = tmp[fid]
                matrix[tmp[iid]] = foo
        fpkm_f.close()

    #OUTPUT:
    print(",".join([options.name]+samples))
    for k in matrix.keys():
        print(",".join([k]+matrix[k]))


if __name__=='__main__':
    main()


