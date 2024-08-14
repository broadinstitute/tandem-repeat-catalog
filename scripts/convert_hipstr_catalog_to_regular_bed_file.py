import argparse
import gzip
import re
import os

p = argparse.ArgumentParser()
p.add_argument("-o", "--output-path", help="File path of the uncrompressed output bed file")
p.add_argument("hipstr_catalog_bed", help="Path of .bed file")
args = p.parse_args()

print(f"Processing HipSTR catalog: {args.hipstr_catalog_bed}")

output_path = args.output_path
if not output_path:
    output_path = re.sub(".bed(.gz)?$", "", args.hipstr_catalog_bed) + ".adjusted.bed"
elif output_path.endswith(".gz"):
    output_path = output_path.replace(".gz", "")

fopen = gzip.open if args.hipstr_catalog_bed.endswith("gz") else open
with fopen(args.hipstr_catalog_bed, "rt") as f, open(output_path, "wt") as fo:
    total_counter = 0
    for i, line in enumerate(f):
        total_counter += 1
        fields = line.strip("\n").split("\t")

        # convert start to 0-based coords standard for bed files
        chrom = fields[0]
        start_0based = int(fields[1]) - 1
        end = int(fields[2])

        # swap columns 3 and 4 so that it's motif first and then motif size. Drop columns after that.
        motifs = fields[6].split("/") # some records in the HipSTR catalog have motif specified as "CCTG/CCCT/CCTT"
        motif = motifs[0]  # the adjusted catalog will only be used for overlap checks, so only motif length matters

        if not re.match("^[ACGTN]+$", motif):
            raise ValueError(f"ERROR: Motif on line {i+1} contains unexpected chars: {motif}")

        motif_size = len(motif)

        fo.write("\t".join(map(str, [chrom, start_0based, end, motif, motif_size])) + "\n")

print(f"Parsed {total_counter:,d} loci from {args.hipstr_catalog_bed}")

os.system(f"bgzip -f {output_path}")
#os.system(f"tabix {output_path}.gz")
print(f"Wrote {total_counter:,d} rows to {output_path}.gz")
