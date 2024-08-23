import argparse
import pandas as pd
from str_analysis.utils.misc_utils import parse_interval
from str_analysis.utils.eh_catalog_utils import parse_motifs_from_locus_structure, get_variant_catalog_iterator

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--output-file", help="output TSV path")
parser.add_argument("expansion_hunter_catalog", help="ExpansionHunter variant catalog in JSON format")
args = parser.parse_args()

if not args.output_file:
	args.output_file = args.expansion_hunter_catalog.replace(".json", ".locus_id_map.tsv")

output_rows = []
for i, record in enumerate(get_variant_catalog_iterator(args.expansion_hunter_catalog)):
	motifs = parse_motifs_from_locus_structure(record["LocusStructure"])
	if isinstance(record["ReferenceRegion"], list):
		reference_regions = record["ReferenceRegion"]
		locus_ids = record["VariantId"]
	else:
		reference_regions = [record["ReferenceRegion"]]
		locus_ids = [record["LocusId"]]

	if len(motifs) != len(reference_regions):
		raise ValueError(f"LocusStructure motif count != ReferenceRegion list length in "
						 f"variant catalog record #{i+1}: {record}")
	if len(motifs) != len(locus_ids):
		raise ValueError(f"LocusStructure motif count != the LocusId or VariantId count in "
						 f"variant catalog record #{i+1}: {record}")

	for motif, reference_region, locus_id in zip(motifs, reference_regions, locus_ids):
		chrom, locus_start_0based, locus_end_1based = parse_interval(reference_region)
		chrom = chrom.replace("chr", "")
		output_rows.append({
			"Name": locus_id,
			"ID": f"{chrom}-{locus_start_0based}-{locus_end_1based}-{motif}",
		})

df = pd.DataFrame(output_rows)
df.to_csv(args.output_file, sep="\t", index=False)
print(f"Wrote {len(df)} rows to {args.output_file}")



