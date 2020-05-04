#!/usr/bin/env python
"""Script to collect the virusseq results across all samples. outputs to stdout
SampleID, TranscriptID, Counts, FPKM
"""

import os
import sys
from optparse import OptionParser

def main():
    usage = "USAGE: %prog -f [FPKM FILE_1] -f [FPKM FILE_2] ...-f [FPKM FILE_N]-c [LIST of read counts1] -c [LIST of read counts2] ... -c [LIST of read countsN]\n\n **ASSUMES that the FPKM and READ counts correspond**"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-f", "--fpkms", action="append", help="list of virusseq.filtered.gtf files")
    optparser.add_option("-c", "--counts", action="append", help="list of virus.ReadsPerGene.out.tab")
    (options, args) = optparser.parse_args(sys.argv)

    if not options.fpkms or not options.counts:
        optparser.print_help()
        sys.exit(-1)

    if len(options.fpkms) != len(options.counts):
        print("discrepant number of fpkm and counts files")
        sys.exit(-1)

    #TRY to infer the SAMPLE NAMES--SAMPLE.virusseq.ReadsPerGene.out.tab
    sampleIDs=[n.strip().split("/")[-1].split('.')[0] for n in options.counts]

    print(",".join(["SampleID","TranscriptID","Counts","FPKM"]))

    for (i, fpkm_f) in enumerate(options.fpkms):
        #READ in fpkm
        f = open(fpkm_f)
        _dict = {}
        for l in f:
            #GET TRANSCRIPT ID (elm 0) and FPKM (elm 1)
            tmp = [tuple(e.strip().split()) for e in l.strip().split(";")]
            transcript_id = tmp[0][1]
            fpkm = tmp[1][1]
            #EVALS will get rid of the "quotes"
            _dict[eval(transcript_id)] = [0, eval(fpkm)]
        #print(_dict)
        f.close()

        #READ in counts
        #NOTE: must be in the FPKM file!!!
        f = open(options.counts[i])
        for l in f:
            tmp = l.strip().split('\t')
            if (tmp[0] in _dict) and int(tmp[1]) > 0:
                _dict[tmp[0]][0] += int(tmp[1])

        #COMPOSE TABLE
        table = sorted(_dict.items(), key=lambda x: x[1][0], reverse=True)

        #OUTPUT:
        for r in table:
            print(",".join([sampleIDs[i], r[0], str(r[1][0]), str(r[1][1])]))

if __name__=='__main__':
    main()


