import argparse
import gzip
import ijson
import simplejson as json
import tqdm


def main():
	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("--verbose", action="store_true")
	parser.add_argument("--show-progress-bar", action="store_true", help="Show a progress bar")
	parser.add_argument("--output-json-path", help="Path of the output JSON file")
	parser.add_argument("variation_clusters_bed_path", help="Path of the variation clusters BED file")
	parser.add_argument("catalog_json_path", help="Path of the JSON catalog to annotate")
	args = parser.parse_args()

	if not args.output_json_path:
		args.output_json_path = args.catalog_json_path.replace(".json", ".with_variation_clusters.json")

	# parse the variation clusters bed file
	#locus_id_to_variation_cluster_id = {}
	locus_id_to_variation_cluster_interval = {}

	if args.verbose:
		print(f"Parsing {args.variation_clusters_bed_path}")
	fopen = gzip.open if args.variation_clusters_bed_path.endswith("gz") else open
	with fopen(args.variation_clusters_bed_path, "rt") as f:
		if args.show_progress_bar:
			f = tqdm.tqdm(f, unit=" records", unit_scale=True)
		for line in f:
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
			for locus_id in info_fields_dict["ID"].split(","):
				locus_id_to_variation_cluster_interval[locus_id] = f"{chrom}:{start_0based}-{end_1based}"
				#locus_id_to_variation_cluster_id[locus_id] = info_fields_dict["ID"]

	total_variation_clusters = len(locus_id_to_variation_cluster_interval)
	if args.verbose:
		print(f"Parsed {len(locus_id_to_variation_cluster_interval):,d} variation clusters from {args.variation_clusters_bed_path}")

	fopen = gzip.open if args.catalog_json_path.endswith("gz") else open
	with fopen(args.catalog_json_path, "rt") as f:
		f2open = gzip.open if args.output_json_path.endswith("gz") else open
		with f2open(args.output_json_path, "wt") as f2:
			input_counter = skipped_counter = output_counter = 0
			iterator = ijson.items(f, "item")
			if args.show_progress_bar:
				iterator = tqdm.tqdm(iterator, unit=" records", unit_scale=True)
			f2.write("[")
			for i, record in enumerate(iterator):
				locus_id = record["LocusId"]
				input_counter += 1
				if locus_id not in locus_id_to_variation_cluster_interval:
					skipped_counter += 1
					#if args.verbose:
					#	print(f"WARNING: locus_id {locus_id} not found in variation cluster catalog")
					continue

				#record["VariationClusterId"] = locus_id_to_variation_cluster_id[locus_id]
				record["VariationCluster"] = locus_id_to_variation_cluster_interval[locus_id]
				#del locus_id_to_variation_cluster_id[locus_id]
				del locus_id_to_variation_cluster_interval[locus_id]

				output_counter += 1
				if i > 0:
					f2.write(",")
				f2.write(json.dumps(record, f2, use_decimal=True, indent=4))
				f2.write("\n")

			f2.write("]")

	if args.verbose:
		print(f"Skipped {skipped_counter:,d} out of {input_counter:,d} ({skipped_counter/input_counter:.1%}) records "
			  f"from {args.catalog_json_path}")
		print(f"Used {total_variation_clusters - len(locus_id_to_variation_cluster_interval):,d} out of {total_variation_clusters:,d} "
			  f"({(total_variation_clusters - len(locus_id_to_variation_cluster_interval))/total_variation_clusters:.1%}) of "
			  f"variation clusters from {args.variation_clusters_bed_path}")
		print(f"Wrote {output_counter:,d} out of {input_counter:,d} ({output_counter/input_counter:.1%}) records "
			  f"to {args.output_json_path}")


if __name__ == "__main__":
	main()

