"""Convert the simple repeat track from UCSC to BED catalog format for comparison.

Source:
https://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=rep&hgta_track=simpleRepeat&hgta_table=simpleRepeat&hgta_doSchema=describe+table+schema
"""

import argparse
import gzip
import json
import os
import tqdm

def main():
	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("--verbose", action="store_true")
	parser.add_argument("-o", "--output-bed", help="Output BED file path")
	parser.add_argument("--show-progress-bar", action="store_true", help="Show a progress bar")
	parser.add_argument("simple_repeat_track_txt")
	args = parser.parse_args()

	if not os.path.isfile(args.simple_repeat_track_txt):
		parser.error(f"File not found: {args.simple_repeat_track_txt}")

	print(f"Converting {args.simple_repeat_track_txt} to BED format")
	if not args.output_bed:
		args.output_bed = args.simple_repeat_track_txt.replace(".txt", ".bed")
	if args.output_bed.endswith("gz"):
		args.output_bed = args.output_bed.replace(".gz", "")

	counter = 0
	supercontig_loci_counter = 0

	fopen = gzip.open if args.simple_repeat_track_txt.endswith("gz") else open
	with fopen(args.simple_repeat_track_txt, "rt") as f:
		if args.show_progress_bar:
			f = tqdm.tqdm(f, unit=" records", unit_scale=True)

		with open(args.output_bed, "wt") as out:
			for line in f:
				counter += 1
				fields = line.strip("\n").split("\t")
				chrom = fields[1]
				if "_" in chrom:
					supercontig_loci_counter += 1
					continue
				start_0based = int(fields[2])
				end = int(fields[3])
				motif = fields[-1]
				out.write(f"{chrom}\t{start_0based}\t{end}\t{motif}\n")

	os.system(f"bgzip -f {args.output_bed}")
	if supercontig_loci_counter:
		print(f"Skipped {supercontig_loci_counter:,d} out of {counter:,d} loci "
			  f"({supercontig_loci_counter/counter:.2%}) because they are on supercontigs")
	print(f"Wrote {counter:,d} rows to {args.output_bed}.gz")

if __name__ == "__main__":
	main()