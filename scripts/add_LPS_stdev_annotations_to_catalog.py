"""Add longest pure segment (LPS) annotations to catalog:

LPSLengthStdevFromHPRC100 is the standard deviation of the length of the longest pure segment (LPS) detected at
	a tandem repeat locus in 100 high-coverage long-read (LR) samples from the HPRC  (Range: 0 to 1226.0)
LPSMotifFractionFromHPRC100 is the fraction of alleles at a tandem repeat locus in which the given motif
	composed the longest pure segment (LPS) among 100 high-coverage long-read (LR) samples from the HPRC.
"""

import argparse
import collections
import gzip
import ijson
import os
import pandas as pd
import simplejson as json
import sys
import tqdm

from str_analysis.utils.misc_utils import parse_interval
from str_analysis.utils.eh_catalog_utils import parse_motifs_from_locus_structure
"""
Expected columns in lps table:
'TRID', 'longestPureSegmentMotif', 'N_motif', '0thPercentile',
       '1stPercentile', '5thPercentile', '10thPercentile', '15thPercentile',
       '20thPercentile', '25thPercentile', '30thPercentile', '35thPercentile',
       '40thPercentile', '45thPercentile', '50thPercentile', '55thPercentile',
       '60thPercentile', '65thPercentile', '70thPercentile', '75thPercentile',
       '80thPercentile', '85thPercentile', '90thPercentile', '95thPercentile',
       '99thPercentile', '99.9thPercentile', '100thPercentile', 'MAD', 'Mean',
       'Stdev'
"""

def main():
	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("--known-pathogenic-loci-json-path", required=True, help="Path of ExpansionHunter catalog "
						"containing known pathogenic loci. This is used to retrieve the original locus boundaries for "
						"these loci since their IDs don't contain these coordinates the way that IDs of other loci do.")
	parser.add_argument("--show-progress-bar", action="store_true", help="Show a progress bar")
	parser.add_argument("--output-catalog-json-path",
						help="Path of the output catalog JSON file that includes variation cluster annotations")
	parser.add_argument("lps_table", help="Path of the LPS data table", default="HPRC_100_LongestPureSegmentQuantiles.txt.gz")
	parser.add_argument("catalog_json_path", help="Path of the JSON catalog to annotate")
	args = parser.parse_args()

	for path in args.lps_table, args.catalog_json_path, args.known_pathogenic_loci_json_path:
		if not os.path.isfile(path):
			parser.error(f"{path} file not found")

	if not args.output_catalog_json_path:
		args.output_catalog_json_path = args.catalog_json_path.replace(".json", ".with_LPS_annotations.json")

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

	print(f"Parsed {len(known_pathogenic_reference_regions_lookup)} known pathogenic loci")
	print(f"Parsing {args.lps_table}")
	df = pd.read_table(args.lps_table)
	missing_columns = {"TRID", "longestPureSegmentMotif", "N_motif", "Stdev"} - set(df.columns)
	if missing_columns:
		parser.error(f"{args.lps_table} is missing expected columns: {missing_columns}")

	before = len(df)
	df = df[~df["longestPureSegmentMotif"].isna() & ~df["Stdev"].isna() & ~df["N_motif"].isna()]
	print(f"Filtered out {before - len(df):,d} out of {before:,d} ({(before - len(df)) / before:.1%}) records with missing values")
	
	# sum the N_motif column across all rows with the same TRID
	TRID_to_N_motif_sum_lookup = dict(df.groupby("TRID")["N_motif"].sum())

	annotation_lookup = {}
	row_iterator = df.iterrows()
	if args.show_progress_bar:
		row_iterator = tqdm.tqdm(row_iterator, total=len(df), unit=" records", unit_scale=True)
		
	for _, row in row_iterator:
		TRID = row["TRID"]
		motif_fraction_string = f"{row['longestPureSegmentMotif']}: {row['N_motif']}/{TRID_to_N_motif_sum_lookup[TRID]}"
		# convert stdev in bp to stdev in repeat units
		stdev_string = "%0.3f" % float(row["Stdev"] / len(row['longestPureSegmentMotif']))

		for locus_id in TRID.split(","):
			if locus_id in known_pathogenic_reference_regions_lookup:
				reference_region, motif = known_pathogenic_reference_regions_lookup[locus_id]
			else:
				assert locus_id.count("-") == 3
				chrom, start_0based, end, motif = locus_id.split("-")
				reference_region = f"{chrom}:{start_0based}-{end}"
			
			if motif != row["longestPureSegmentMotif"]:
				continue
				
			annotation_lookup[locus_id] = {
				"LPSLengthStdevFromHPRC100": stdev_string,
				"LPSMotifFractionFromHPRC100": motif_fraction_string,
			}

	input_locus_counter = annotated_locus_counter = 0
	print(f"Adding LPS annotations to {args.catalog_json_path}")
	fopen = gzip.open if args.catalog_json_path.endswith("gz") else open
	with fopen(args.catalog_json_path, "rt") as f:
		f2open = gzip.open if args.output_catalog_json_path.endswith("gz") else open
		with f2open(args.output_catalog_json_path, "wt") as f2:
			iterator = ijson.items(f, "item")
			if args.show_progress_bar:
				iterator = tqdm.tqdm(iterator, unit=" records", unit_scale=True)
			f2.write("[")
			for i, record in enumerate(iterator):
				locus_id = record["LocusId"]
				input_locus_counter += 1
				if locus_id in annotation_lookup:
					record.update(annotation_lookup[locus_id])
					annotated_locus_counter += 1
				if i > 0:
					f2.write(", ")
				f2.write(json.dumps(record, f2, use_decimal=True, indent=4))
			f2.write("]")

	print(f"Annotated {annotated_locus_counter:,d} out of {input_locus_counter:,d} loci")
	print(f"Wrote annotated catalog to {args.output_catalog_json_path}")

if __name__ == "__main__":
	main()

