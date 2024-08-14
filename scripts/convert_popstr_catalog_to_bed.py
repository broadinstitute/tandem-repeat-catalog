"""This script takes the popSTR chrXXXmarkerInfo.gz files that are available at
https://github.com/DecodeGenetics/popSTR and converts/combines them into a BED file.

Example row from the chrXXXmarkerInfo.gz input file format:
$1                     chrom : chr22
$2           startCoordinate : 10517060
$3             endCoordinate : 10517070
$4               repeatMotif : ATGAG
$5         numOfRepeatsInRef : 2.2
$6   1000refBasesBeforeStart : AAAAGTCCATCATCAAATGAACAGATGAAGAAAATATGGTATATGTGTGTGTGTGTATATATATATATGTATGGTATATATATGTATGGTGTATATATATATATATGTATGATATATATAGTATGGCATATATATGTATGGTGTATATATATGGTATATTTATGTGTATATATATGTATATATGTATATATATGTATATATACATACACACACAGAATGGAATATTAGTCAGCCTTCAAAAGGAAAATTCTGTCGTATTTCAACATGTATCAAGCTTAAGGATATTGTGCTAAGTGAAATAAGCCAGACACAAAGACAAATATATCGTGATTCCATTTATATGATGTATTGAAAGTAGCCAAACACATGGAAACACAAGATAAAATGGTAGTGGTCAGGGCCTGGAGGAAACAGGAAATCTGGAGTTGCTGTTCACCAGGTGTAGAGTTTCAGTCAAGCAAGATAAAAACATTCTAGATTTCTGCTGTACAACAATGTGTATATCATTAACAAAATGTTCTGCAAACTTAAAATTTTGTTAAGATGGTAGATTTTTTTGTTATGTGTTTTTTAATTACAAAAATTTCTGTCTGTATTTAGTTTACATTTTAGTAAGAAAAGACAAATAACCTAATAAATGGTAATGATAAATGCTGTAAGGATATCTAGAGCAATAAAAGAAGTTATGGATATGGGAGTATAATTTTAGATAGTGAAGATTTCTGTATTCAAATGCCACATGCAAAAAGGACTAAGGGGAATGAGGAGATGAGTCATATGGAATGCCCAAGACACAGAAAAAGGCAGGCAGAGAAAATAACAAGGATAAAGACACTGAAGTAAAATCATGCTTTCTATGTTTCAAAACAGCAGCAAGTACACTACTGTGATAGAGCAGGAGTGACCAATGGGAGGCTGGAGATTTTGTCAGAGATATTGTCAAGGCTCAGGATCGTACAGGGACTTGTAAGCCTGGAAAGCACTTTGAATTTTATTCAGA
$7      1000refBasesAfterEnd : AGCCATTGAAAAGTTTTTAAGCAGATGAGTAAAATAATCCACCTTGTATTTTAAGAGGAGCATTCTACCTTCTCTGTGGAATAGAGAGGTGGAAGGGAAAGCTTGAAGCAGAGAGAGCAGTGAAGAGTGTGCTGTAATATTCTTATGTGAGAAATAGTGGAGGGAATGAGAGGTGGTCAGCCTTAAACTGCCATTTGCTCTCTGTATCAGGGCTCAGGGACTTTCAGACTCTCCAGGGATTCCTTACAGTTTTTACATTTGTTTCTCAGTTGCAAACATAATCTCTTCACTCTATGAGAACTCTAGATCTGAATCCTTGTTATGAGTCAGGAGTCCTCTCTAGTTTACTCTCTTGTCAGTAACTAGACTTGAAATGTTTTGATATTAATTGATATAGTAATAAGATTGGTTAGAGAAATAGCAAAGAGAGAGCATCCCCATCCTATGACCATATCAGCACCAGAAGAGAAAAACACATCTACACAGTTTTTCCCTTGGCATAGGTCCTGGTATTCTGTTAGGGAACAGGTTATGAAGGAAATACAAAGGGTCTTTGTACTTACTTTCAGAATGTATTTTCTTTAACATGAAAAGAATCCAAGGCCTTTTTGCTTCTAATTGCTTTTTGTGTATCTACTACCATCCCTTGCTAAATTATTGATAGGTTTCCTCAAATCTCGGCATGATGTCCTACATTCTAAATTTTCAATAGCTGAAAATTTCACCTTTTCAGTGCCTTCAAGTTTATCTCAGTAAAAAGTTGAGAAAGACTGTAATAGAGTTATTTAATCAGATTTTTTTCATCTACCATAATTTTTGAATAAGGAAAAACAGCAATACTTTTTCTCCTTACTTGGCAAGTAATTTTCATAGAGAGGAAAAAAACAATCAAAACAGGTACAAAATGTAACAAAACCAAAGGACCATGTGAGGTGAAATTTAAAATGAGAAAAATGTCCACAGTACTTTGGGCAATGCAACTCCTGAGAAATAGTAACTC
$8          repeatSeqFromRef : ATGAGATGAGA
$9              minFlankLeft : 4
$10            minFlankRight : 4
$11             repeatPurity : 1.00
$12          defaultSlippage : 0.017275
$13           defaultStutter : 0.946913
$14         fractionAinMotif : 0.4
$15         fractionCinMotif : 0
$16         fractionGinMotif : 0.4
$17         fractionTinMotif : 0.2
"""
import argparse
import gzip
import os

parser = argparse.ArgumentParser()
parser.add_argument("--output-bed", default="popstr_catalog_v2.bed", help="Output bed file path")
parser.add_argument("popstr_marker_info_files_gz", nargs="+",
                    help="One or more of the chrXXXmarkerInfo.gz files downloaded/cloned from "
                         "https://github.com/DecodeGenetics/popSTR ")
args = parser.parse_args()

if not args.output_bed.endswith(".bed"):
    parser.error(f"Output file path must have a .bed suffix: {args.output_bed}")

# parse input files
output_rows = []
for path in args.popstr_marker_info_files_gz:
    print(f"Parsing {path}")
    with gzip.open(path, "rt") as input_file:
        for line in input_file:
            fields = line.strip().split()
            chrom = fields[0]
            start_1based = int(fields[1])
            end = int(fields[2])
            motif = fields[3]
            output_rows.append((chrom, start_1based - 1, end, motif, "."))

print(f"Parsed {len(output_rows):,d} rows from {len(args.popstr_marker_info_files_gz)} input file(s)")

# write all rows to output bed
with open(args.output_bed, "w") as f:
    for row in sorted(output_rows):
        f.write("\t".join(map(str, row)) + "\n")

print(f"Wrote {len(output_rows):,d} rows to {args.output_bed}")

# compress and index output bed
def run(command):
    print(command)
    os.system(command)

run(f"bgzip -f {args.output_bed}")
run(f"tabix -f -p bed {args.output_bed}.gz")

print("Done")