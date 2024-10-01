
import argparse
import collections
import gzip
import os
import pandas as pd
import simplejson as json
import sys
import tqdm

from str_analysis.utils.eh_catalog_utils import parse_motifs_from_locus_structure
from str_analysis.utils.find_repeat_unit import find_repeat_unit_without_allowing_interruptions
from str_analysis.utils.misc_utils import parse_interval


stats = collections.Counter()

def fix_variation_cluster_id(variation_cluster_id, known_pathogenic_reference_regions_lookup):
	new_ids = []
	updated = False
	stats["total"] += 1
	for locus_id in variation_cluster_id.split(","):
		if locus_id in known_pathogenic_reference_regions_lookup:
			reference_region, motif = known_pathogenic_reference_regions_lookup[locus_id]
			chrom, start_0based, end_1based = parse_interval(reference_region)
			chrom = chrom.replace("chr", "")
			new_ids.append(f"{chrom}-{start_0based}-{end_1based}-{motif}")
			updated = True
		else:
			assert locus_id.count("-") == 3, f"Unexpected locus ID format: {locus_id}"
			chrom, start_0based, end_1based, motif = locus_id.split("-")

			simplified_motif, num_repeats, _ = find_repeat_unit_without_allowing_interruptions(
				motif, allow_partial_repeats=False)

			chrom = chrom.replace("chr", "")
			new_ids.append(f"{chrom}-{start_0based}-{end_1based}-{simplified_motif}")
			if simplified_motif != motif:
				#print(f"Updating {locus_id} to {new_ids[-1]}")
				updated=True

		if updated:
			stats["locus_id_updated"] += 1

	if updated:
		stats["variation_cluster_id_updated"] += 1

	return ",".join(new_ids)

def main():
	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("--known-pathogenic-loci-json-path", required=True, help="Path of ExpansionHunter catalog "
						"containing known pathogenic loci. This is used to retrieve the original locus boundaries for "
						"these loci since their IDs don't contain these coordinates the way that IDs of other loci do.")
	parser.add_argument("--verbose", action="store_true")
	parser.add_argument("--show-progress-bar", action="store_true", help="Show a progress bar")
	parser.add_argument("input_path", help="Path of the variation clusters BED file or the LPS data table")
	args = parser.parse_args()

	for path in args.input_path, args.known_pathogenic_loci_json_path:
		if not os.path.isfile(path):
			parser.error(f"File not found: {path}")

	print(f"Parsing {args.known_pathogenic_loci_json_path}")
	fopen = gzip.open if args.known_pathogenic_loci_json_path.endswith("gz") else open
	with fopen(args.known_pathogenic_loci_json_path, "rt") as f:
		known_pathogenic_loci = json.load(f)
		known_pathogenic_reference_regions_lookup = {}
		for locus in known_pathogenic_loci:
			motifs = parse_motifs_from_locus_structure(locus["LocusStructure"])
			if isinstance(locus["ReferenceRegion"], list):
				assert isinstance(locus["VariantId"], list)
				assert len(locus["ReferenceRegion"]) == len(locus["VariantId"])
				assert len(locus["ReferenceRegion"]) == len(motifs)
				for variant_id, reference_region, motif in zip(locus["VariantId"], locus["ReferenceRegion"], motifs):
					known_pathogenic_reference_regions_lookup[variant_id] = (reference_region, motif)
			else:
				known_pathogenic_reference_regions_lookup[locus["LocusId"]] = (locus["ReferenceRegion"], motifs[0])

	if args.input_path.endswith(".bed") or args.input_path.endswith(".bed.gz"):
		print(f"Parsing BED file: {args.input_path}")

		output_path = args.input_path.replace(".bed", ".fixed_variation_cluster_ids.bed")
		assert output_path != args.input_path, f"Unexpected input path: {args.input_path}"
		print(f"Writing to {output_path}")
		fopen = gzip.open if args.input_path.endswith("gz") else open
		with fopen(args.input_path, "rt") as f, fopen(output_path, "wt") as out_f:
			counter = 0
			for line in f:
				fields = line.strip().split("\t")
				assert len(fields) >= 4, f"Expected at least 4 fields in BED file, got {len(fields)}: {line}"
				info_fields = fields[3]
				info_fields_dict = {}
				for key_value in info_fields.split(";"):
					key_value = key_value.split("=")
					if len(key_value) != 2:
						print(f"WARNING: skipping invalid key-value pair '{key_value}' in line {fields}")
						continue
					key, value = key_value
					info_fields_dict[key] = value

				new_variation_cluster_id = fix_variation_cluster_id(
					info_fields_dict["ID"], known_pathogenic_reference_regions_lookup)
				if new_variation_cluster_id != info_fields_dict["ID"]:
					info_fields_dict["ID"] = new_variation_cluster_id
					fields[3] = ";".join(f"{k}={v}" for k, v in info_fields_dict.items())

				counter += 1
				out_f.write("\t".join(fields) + "\n")
			print(f"Wrote {counter:,d} lines to {output_path}")
	elif (
		args.input_path.endswith(".txt") or args.input_path.endswith(".txt.gz")
		or args.input_path.endswith(".tsv") or args.input_path.endswith(".tsv.gz")
	):
		if args.input_path.endswith(".txt") or args.input_path.endswith(".txt.gz"):
			output_path = args.input_path.replace(".txt", ".fixed_variation_cluster_ids.txt")
		elif args.input_path.endswith(".tsv") or args.input_path.endswith(".tsv.gz"):
			output_path = args.input_path.replace(".tsv", ".fixed_variation_cluster_ids.tsv")
		else:
			assert False, f"Unexpected file format: {args.input_path}"
			
		assert output_path != args.input_path, f"Unexpected input path: {args.input_path}"
		print(f"Parsing table: {args.input_path}")
		df = pd.read_table(args.input_path)
		df["TRID"] = df["TRID"].apply(lambda trid: fix_variation_cluster_id(
					trid, known_pathogenic_reference_regions_lookup))
		df.to_csv(output_path, sep="\t", index=False)
		print(f"Wrote {len(df):,d} rows to {output_path}")
	else:
		parser.error(f"Unsupported file format: {args.input_path}")

	if stats:
		print(f"Updated {stats['variation_cluster_id_updated']:,d} out of {stats['total']:,d} "
			  f"({stats['variation_cluster_id_updated']/stats['total']:.2%}) variation cluster IDs")


if __name__ == "__main__":
	main()