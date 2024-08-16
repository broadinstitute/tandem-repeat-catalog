import argparse
import collections
import gzip
import ijson
import os
import simplejson as json
import sys
import tqdm

from str_analysis.utils.misc_utils import parse_interval

ALMOST_NO_CHANGE_TO_BOUNDARIES_THRESHOLD = 1

def main():
	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("--verbose", action="store_true")
	parser.add_argument("--known-pathogenic-loci-json-path", required=True, help="Path of ExpansionHunter catalog "
						"containing known pathogenic loci. This is used to retrieve the original locus boundaries for "
						"these loci since their IDs don't contain these coordinates the way that IDs of other loci do.")
	parser.add_argument("--show-progress-bar", action="store_true", help="Show a progress bar")
	parser.add_argument("--output-catalog-json-path", 
						help="Path of the output catalog JSON file that includes variation cluster annotations")
	parser.add_argument("--output-variation-clusters-bed-path", 
						help="Path of the output variation clusters BED file that excludes variation clusters whose "
						"boundaries exactly match those of loci found in the input catalog")
	parser.add_argument("variation_clusters_bed_path", help="Path of the variation clusters BED file")
	parser.add_argument("catalog_json_path", help="Path of the JSON catalog to annotate")
	args = parser.parse_args()

	if not args.output_catalog_json_path:
		args.output_catalog_json_path = args.catalog_json_path.replace(".json", ".with_variation_clusters.json")
	if not args.output_variation_clusters_bed_path:
		args.output_variation_clusters_bed_path = args.variation_clusters_bed_path.replace(".bed", ".excluding_duplicates_of_catalog_loci.bed")
	if not os.path.isfile(args.known_pathogenic_loci_json_path):
		parser.error(f"File not found: {args.known_pathogenic_loci_json_path}")

	with open(args.known_pathogenic_loci_json_path, "rt") as f:
		known_pathogenic_loci = json.load(f)
		known_pathogenic_reference_regions_lookup = {}
		for locus in known_pathogenic_loci:
			if isinstance(locus["ReferenceRegion"], list):
				assert isinstance(locus["VariantId"], list)
				assert len(locus["ReferenceRegion"]) == len(locus["VariantId"])
				for variant_id, reference_region in zip(locus["VariantId"], locus["ReferenceRegion"]):
					known_pathogenic_reference_regions_lookup[variant_id] = reference_region
			else:
				known_pathogenic_reference_regions_lookup[locus["LocusId"]] = locus["ReferenceRegion"]

	locus_id_to_variation_cluster_interval = {}
	locus_id_to_variation_cluster_size_difference_from_simple_repeat_boundaries = {}
	size_diff_histogram = collections.Counter()
	input_variation_clusters_counter = 0
	input_locus_ids_counter = 0
	almost_no_change_to_boundaries = 0
	output_variation_clusters_counter = 0
	if args.verbose:
		print(f"Parsing {args.variation_clusters_bed_path}")

	fopen = gzip.open if args.variation_clusters_bed_path.endswith("gz") else open
	fopen2 = gzip.open if args.output_variation_clusters_bed_path.endswith("gz") else open
	with fopen(args.variation_clusters_bed_path, "rt") as f, fopen2(args.output_variation_clusters_bed_path, "wt") as output_bed_file:
		if args.show_progress_bar:
			f = tqdm.tqdm(f, unit=" records", unit_scale=True)

		for line in f:
			input_variation_clusters_counter += 1
			fields = line.strip("\n").split("\t")
			chrom = fields[0]
			start_0based = int(fields[1])
			end_1based = int(fields[2])
			info_fields = fields[3]

			info_fields_dict = {}
			for key_value in info_fields.split(";"):
				key_value = key_value.split("=")
				if len(key_value) != 2:
					print(f"WARNING: skipping invalid key-value pair '{key_value}' in line {fields}")
					continue
				key, value = key_value
				info_fields_dict[key] = value

			variation_cluster_differs_from_simple_repeat = False
			region = f"{chrom.replace('chr', '')}:{start_0based}-{end_1based}"
			for locus_id in info_fields_dict["ID"].split(","):
				input_locus_ids_counter += 1
				if locus_id in known_pathogenic_reference_regions_lookup:
					region2 = known_pathogenic_reference_regions_lookup[locus_id]
					original_chrom, original_start_0based, original_end_1based = parse_interval(region2)
				elif locus_id.count("-") == 3:
					original_chrom, original_start_0based, original_end_1based, _ = locus_id.split("-")
					region2 = f"{original_chrom}:{original_start_0based}-{original_end_1based}"
				else:
					raise ValueError(f"Unexpected locus_id '{locus_id}'")

				original_start_0based = int(original_start_0based)
				original_end_1based = int(original_end_1based)

				size_diff = (end_1based - start_0based) - (original_end_1based - original_start_0based)
				if size_diff < 0:
					print(f"WARNING: {locus_id} variation cluster {region} is smaller than the original locus {region2}")

				if abs(size_diff) <= ALMOST_NO_CHANGE_TO_BOUNDARIES_THRESHOLD:
					almost_no_change_to_boundaries += 1
				else:
					variation_cluster_differs_from_simple_repeat = True
					locus_id_to_variation_cluster_interval[locus_id] = region
					locus_id_to_variation_cluster_size_difference_from_simple_repeat_boundaries[locus_id] = size_diff
					size_diff_histogram[size_diff] += 1

			if variation_cluster_differs_from_simple_repeat:
				output_bed_file.write(line)
				output_variation_clusters_counter += 1

	total_locus_ids = len(locus_id_to_variation_cluster_interval)
	if args.verbose:
		print(f"Found {total_locus_ids:,d} out of {input_locus_ids_counter:,d} "
			  f"({total_locus_ids/input_locus_ids_counter:.1%}) locus IDs in the input catalog "
			  f"were in variation clusters")
		print(f"Found {almost_no_change_to_boundaries:,d} out of {input_variation_clusters_counter:,d} "
			  f"({almost_no_change_to_boundaries/input_variation_clusters_counter:.1%}) "
			  f"variation clusters did not change the original locus boundaries "
			  f"by more than {ALMOST_NO_CHANGE_TO_BOUNDARIES_THRESHOLD} bases")
		print(f"Wrote {output_variation_clusters_counter:,d} out of {input_variation_clusters_counter:,d} "
			  f"({output_variation_clusters_counter/input_variation_clusters_counter:.1%}) variation clusters "
			  f"to {args.output_variation_clusters_bed_path}")
		print(f"Found {output_variation_clusters_counter:,d} out of {input_variation_clusters_counter:,d} "
			  f"({output_variation_clusters_counter/input_variation_clusters_counter:.1%}) "
			  f"variation clusters that differed from locus definitions in the "
			  f"input catalog {args.catalog_json_path}")


	fopen = gzip.open if args.catalog_json_path.endswith("gz") else open
	with fopen(args.catalog_json_path, "rt") as f:
		f2open = gzip.open if args.output_catalog_json_path.endswith("gz") else open
		with f2open(args.output_catalog_json_path, "wt") as f2:
			input_locus_counter = locus_without_variation_cluster_counter = output_locus_counter = 0
			iterator = ijson.items(f, "item")
			if args.show_progress_bar:
				iterator = tqdm.tqdm(iterator, unit=" records", unit_scale=True)
			f2.write("[")
			for i, record in enumerate(iterator):
				locus_id = record["LocusId"]
				input_locus_counter += 1
				if locus_id in locus_id_to_variation_cluster_interval:
					#record["VariationClusterId"] = locus_id_to_variation_cluster_id[locus_id]
					record["VariationCluster"] = locus_id_to_variation_cluster_interval[locus_id]
					record["VariationClusterSizeDiff"] = locus_id_to_variation_cluster_size_difference_from_simple_repeat_boundaries[locus_id]
					#del locus_id_to_variation_cluster_id[locus_id]
					del locus_id_to_variation_cluster_interval[locus_id]
				else:
					locus_without_variation_cluster_counter += 1
					#if args.verbose:
					#	print(f"WARNING: locus_id {locus_id} not found in variation cluster catalog")

				output_locus_counter += 1
				if i > 0:
					f2.write(", ")
				f2.write(json.dumps(record, f2, use_decimal=True, indent=4))
			f2.write("]")

	if args.verbose:
		print(f"{locus_without_variation_cluster_counter:,d} out of {input_locus_counter:,d} "
			  f"({locus_without_variation_cluster_counter/input_locus_counter:.1%}) loci from {args.catalog_json_path} are not in variation clusters")
		print(f"Wrote {output_locus_counter:,d} out of {input_locus_counter:,d} ({output_locus_counter/input_locus_counter:.1%}) records "
			  f"to {args.output_catalog_json_path}")

		#print(size_diff_histogram)

		# plot using seaborn
		import seaborn as sns
		import matplotlib.pyplot as plt
		plt.figure(figsize=(12, 6))
		sns.barplot(x=list(size_diff_histogram.keys()), y=list(size_diff_histogram.values()))
		plt.xlabel("Size difference")
		plt.ylabel("Count")
		plt.title("Size difference between variation clusters and original loci")
		plt.savefig(args.output_variation_clusters_bed_path.replace(".bed.gz", "") + ".size_diff_histogram.png")
		# also create a plot with log y-scale
		plt.yscale("log")
		plt.savefig(args.output_variation_clusters_bed_path.replace(".bed.gz", "") + ".size_diff_histogram.log.png")


if __name__ == "__main__":
	main()

