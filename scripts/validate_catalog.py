"""Run basic sanity checks on the final catalogs prior to release"""

import argparse
import gzip
import json
import ijson
from str_analysis.utils.misc_utils import parse_interval
from tqdm import tqdm

EXPECTED_KEYS_IN_ANNOTATED_CATALOG = {
	"ReferenceRegion": str,
	"LocusStructure": str,
	"VariantType": str,
	"LocusId": str,
	"Source": str,
	"CanonicalMotif": str,
	"GencodeGeneRegion": str,
	#"GencodeGeneName": str,    <== these are optional,
	#"GencodeGeneId": str,
	#"GencodeTranscriptId": str,
	"InterruptionBaseCount": int,
	"FractionPureBases": float,
	"FractionPureRepeats": float,
	"LeftFlankMappability": float,
	"FlanksAndLocusMappability": float,
	"RightFlankMappability": float,
}

def parse_known_pathogenic_loci(known_pathogenic_loci_json_path):
	catalog = {}
	fopen = gzip.open if known_pathogenic_loci_json_path.endswith("gz") else open
	with fopen(known_pathogenic_loci_json_path, "rt") as f:
		known_pathogenic_loci = json.load(f)
		for record in known_pathogenic_loci:
			catalog[record["LocusId"]] = record

	return catalog

def main():
	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("--verbose", action="store_true")
	parser.add_argument("--known-pathogenic-loci-json-path", required=True, help="Path of ExpansionHunter catalog "
						"containing known pathogenic loci. This is used to retrieve the original locus boundaries for "
						"these loci since their IDs don't contain these coordinates the way that IDs of other loci do.")
	parser.add_argument("simple_repeat_catalog_path")
	args = parser.parse_args()

	known_disease_loci_catalog = parse_known_pathogenic_loci(args.known_pathogenic_loci_json_path)

	known_locus_id_to_reference_region = {}
	for locus_id, record in known_disease_loci_catalog.items():
		if isinstance(record["ReferenceRegion"], list):
			for reference_region, variant_id in zip(record["ReferenceRegion"], record["VariantId"]):
				known_locus_id_to_reference_region[variant_id] = reference_region
		else:
			known_locus_id_to_reference_region[locus_id] = record["ReferenceRegion"]

	error_counter = 0

	# confirm that the locus ids and ReferenceRegions are present in the catalog
	locus_ids = set()
	reference_regions = set()
	fopen = gzip.open if args.simple_repeat_catalog_path.endswith("gz") else open
	with fopen(args.simple_repeat_catalog_path, "rt") as f:
		# use ijson to load the catalog iteratively
		for i, record in tqdm(enumerate(ijson.items(f, "item", use_float=True)), unit=" records", unit_scale=True):
			if not record["ReferenceRegion"].startswith("chr"):
				print(f"ERROR: ReferenceRegion {record['ReferenceRegion']} does not start with 'chr'")
				error_counter += 1

			if record["ReferenceRegion"].replace(":", "-").replace("chr", "") not in record["LocusId"] and record["LocusId"] not in known_locus_id_to_reference_region:
				print(f"ERROR: ReferenceRegion {record['ReferenceRegion']} does not match LocusId {record['LocusId']}")
				error_counter += 1

			locus_ids.add(record["LocusId"])
			reference_regions.add(record["ReferenceRegion"])

			for key, value_type in EXPECTED_KEYS_IN_ANNOTATED_CATALOG.items():

				if key not in record or record[key] is None:
					print(f"ERROR: Missing key {key} in record {i}: {record}")
					error_counter += 1
				elif not isinstance(record[key], value_type):
					print(f"ERROR: Expected {key} to be of type {value_type} but got {type(record[key])}")
					error_counter += 1

			if record["GencodeGeneRegion"] != "intergenic":
				for key in ["GencodeGeneName", "GencodeGeneId", "GencodeTranscriptId"]:
					if key not in record or record[key] is None:
						print(f"ERROR: Missing {key} in record {i}: {record}")
						error_counter += 1

	# check that all known disease loci are in the simple repeat catalog
	for locus_id, reference_region in known_locus_id_to_reference_region.items():
		if reference_region not in reference_regions:
			print(f"ERROR: ReferenceRegion {reference_region} not found in the simple repeat catalog")
			error_counter += 1

		if locus_id not in locus_ids:
			print(f"ERROR: {locus_id} not found in the simple repeat catalog")
			error_counter += 1
		else:
			if args.verbose:
				print(f"SUCCESS: {locus_id} found in the simple repeat catalog")


	if error_counter > 0:
		raise ValueError(f"Found {error_counter} errors in {args.simple_repeat_catalog_path}")

if __name__ == "__main__":
	main()