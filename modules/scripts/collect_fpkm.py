"""PART of the immunologymodule
CALLED in generate_cuff_matrix
"""
import sys

gene_list = {}

def path2sample(path):
    # convert "analysis/cufflinks/sampleName/sampleName.genes.fpkm_tracking" into "sampleName"
    parts = path.split("/")
    return parts[-2]
samples = [path2sample(i) for i in sys.argv[1:]]

for fpkm_path in sys.argv[1:]:

    # Remove duplicate gene names
    gene_dedup = set()

    with open(fpkm_path) as f:
        f.readline()

        for line in f:
            cols = line.split("\t")
            # column 4 (index==3) is gene_id, column 10 (index==9)is FPKM
            if cols[3] in gene_dedup:
                continue
            else:
                gene_dedup.add(cols[3])

            if cols[3] in gene_list:
                gene_list[cols[3]].append(cols[9])
            else:
                gene_list[cols[3]] = [cols[9]]

print("GeneName" + "\t" + "\t".join(samples))
for gene in gene_list:
    print(gene, end="\t")
    print("\t".join(gene_list[gene]))
