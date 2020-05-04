IMPORTANT!
How to use these files:
1. replace your genome's original chrM with this chrM, which contains virus
   sequences
2. make a gtf:
   a. remove chrM from your gtf-
      grep -v chrM ../hg38.refseq.gtf > hg38.refseq.noChrM.gtf
   b. add virusSeq.gtf to refseq.noChrM.gtf
3. use the new assembly and new gtf in 1 and 2 to generate a STAR index
