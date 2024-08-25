import argparse
import collections
import gzip
import ijson
import json
import pandas as pd
from str_analysis.utils.file_utils import download_local_copy
from str_analysis.utils.misc_utils import parse_interval
from intervaltree import Interval, IntervalTree

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-o", ]"--output-path", help="Output JSON path for annotated catalog")
parser.add_argument("variant_catalog", help="Variant catalog in JSON format")
args = parser.parse_args()

if not args.output_path:
	args.output_path = args.variant_catalog.replace(".json", ".with_Illumina174k_allele_frequences.json")

def convert_allele_histogram_dict_to_string(allele_histogram_dict):
	data = sorted(allele_histogram_dict.items())
	return ",".join(f"{repeat_number}x:{allele_count}" for repeat_number, allele_count in data)

def get_percentile_from_allele_histogram_dict(allele_histogram_dict, percentile):
	total = sum(allele_histogram_dict.values())
	cutoff = total * percentile
	for repeat_number, count in sorted(allele_histogram_dict.items(), reverse=True):
		total -= count
		if total <= cutoff:
			return repeat_number

def get_stdev_of_allele_histogram_dict(allele_histogram_dict):
	total = sum(allele_histogram_dict.values())
	mean = sum(repeat_number * count for repeat_number, count in allele_histogram_dict.items()) / total
	return (sum((repeat_number - mean) ** 2 * count for repeat_number, count in allele_histogram_dict.items()) / total) ** 0.5


# download illumina table
url = "https://github.com/Illumina/RepeatCatalogs/raw/master/hg38/genotype/1000genomes/1kg.gt.hist.tsv.gz"
print(f"Loading allele frequencies for the Illumina 174k catalog from {url}")
df1 = pd.read_table(download_local_copy(url))
print(f"Parsed {len(df1):,d} rows")
print("Computing histograms")
histograms_from_illumina_174k = {}
stdev_from_illumina_174k = {}
for _, row in df1.iterrows():
	chrom, start_0based, end = row.VariantId.split("_")
	chrom = chrom.replace("chr", "")
	repeat_numbers = [int(x) for x in row.RepeatNumbers.split(",")]
	allele_counts = [int(x) for x in row.AlleleCounts.split(",")]
	if len(repeat_numbers) != len(allele_counts):
		raise ValueError(f"RepeatNumbers and AlleleCounts have different lengths: {row.to_dict()}")

	key = (chrom, int(start_0based), int(end))
	histogram_dict = dict(zip(repeat_numbers, allele_counts))
	histograms_from_illumina_174k[key] = convert_allele_histogram_dict_to_string(histogram_dict)
	stdev_from_illumina_174k[key] = get_stdev_of_allele_histogram_dict(histogram_dict)

print(f"Processed allele frequency histograms for {len(df1):,d} rows")
df1 = None

# download table of genotypes from T2T assemblies
url2 = "gs://str-truth-set-v2/filter_vcf/all_repeats_including_homopolymers/combined/joined_only_pure_repeats.78_samples.variants.tsv.gz"
print(f"Loading allele frequencies for the catalog of polymorphic loci in T2T assemblies from {url2}")
df2 = pd.read_table(download_local_copy(url2))
print(f"Parsed {len(df2):,d} rows")
print("Computing histograms")
allele_columns = [c for c in df2.columns if c.startswith("NumRepeats") and c != "NumRepeatsInReference"]
df2.rename(columns={c: c.replace(":", "_") for c in allele_columns}, inplace=True)
allele_columns = [c.replace(":", "_") for c in allele_columns]

histograms_from_t2t_assemblies = {}
stdev_from_t2t_assemblies = {}
interval_trees_for_t2t_assemblies = collections.defaultdict(IntervalTree)
for row in df2.itertuples():
	chrom, start_1based, end = parse_interval(row.Locus)
	start_0based = start_1based - 1
	chrom = chrom.replace("chr", "")
	key = (chrom, start_0based, end)
	histogram_dict = collections.Counter([
		int(float(getattr(row, c))) for c in allele_columns if not pd.isna(getattr(row, c))
	])
	histograms_from_t2t_assemblies[key] = convert_allele_histogram_dict_to_string(histogram_dict)
	stdev_from_t2t_assemblies[key] = get_stdev_of_allele_histogram_dict(histogram_dict)
	if end > start_0based:
		interval_trees_for_t2t_assemblies[chrom].add(Interval(start_0based, end, data = {
			"HistogramDict": histogram_dict,
			"CanonicalMotif": row.CanonicalMotif,
		}))

print(f"Processed allele frequency histograms for {len(df2):,d} rows")
df2 = None

print(f"Parsing and annotating {args.variant_catalog}")
counters = collections.Counter()
fopen = gzip.open if args.variant_catalog.endswith("gz") else open
fopen2 = gzip.open if args.output_path.endswith("gz") else open
with fopen(args.variant_catalog) as f, fopen2(args.output_path, "wt") as out:
	out.write("[")
	for i, row in enumerate(ijson.items(f, "item", use_float=True)):
		counters["total"] += 1
		if isinstance(row["ReferenceRegion"], list):
			raise ValueError(f"ReferenceRegion is a list in {row.to_dict()}")
		chrom, start_0based, end = parse_interval(row["ReferenceRegion"])
		chrom = chrom.replace("chr", "")
		key = (chrom, start_0based, end)
		if key in histograms_from_illumina_174k:
			counters["found_illumina174_histogram"] += 1
			row["AlleleFrequenciesFromIllumina174k"] = histograms_from_illumina_174k[key]
			row["VariationStdevFromIllumina174k"] = stdev_from_illumina_174k[key]

		if key in histograms_from_t2t_assemblies:
			counters["found_t2t_assemblies_histogram"] += 1
			row["AlleleFrequenciesFromT2TAssemblies"] = histograms_from_t2t_assemblies[key]
			row["VariationStdevFromT2TAssemblies"] = stdev_from_t2t_assemblies[key]
		else:
			# check for overlap with nearby interval
			matching_interval = None
			for interval in interval_trees_for_t2t_assemblies[chrom].overlap(start_0based, end):
				if interval.data["CanonicalMotif"] != row["CanonicalMotif"]:
					continue
				if interval.overlap_size(start_0based, end) < 2*len(row["CanonicalMotif"]):
					continue

				matching_interval = interval
				break

			if matching_interval:
				counters["found_t2t_assemblies_histogram"] += 1
				histogram_dict = matching_interval.data["HistogramDict"]
				locus_boundary_diff = ((matching_interval.end - matching_interval.begin) - (end - start_0based))//len(row["CanonicalMotif"])
				histogram_dict_adjusted = {
					repeat_number - locus_boundary_diff: count for repeat_number, count in histogram_dict.items()
				}
				for repeat_number in histogram_dict_adjusted:
					if repeat_number < 0:
						raise ValueError(f"Negative repeat number {repeat_number} in {histogram_dict_adjusted}")

				row["AlleleFrequenciesFromT2TAssemblies"] = convert_allele_histogram_dict_to_string(histogram_dict_adjusted)
				row["VariationStdevFromT2TAssemblies"] = get_stdev_of_allele_histogram_dict(histogram_dict_adjusted)

				#print(f"Adjusted histogram for {matching_interval} to annotate {row['LocusId']} histogram "
				#	  f"from {convert_allele_histogram_dict_to_string(histogram_dict)} "
				#	  f"to {convert_allele_histogram_dict_to_string(histogram_dict_adjusted)}")


		if i > 0:
			out.write(", ")
		json.dump(row, out, indent=1)
	out.write("]\n")

print(f"Annotated {counters['found_illumina174_histogram']:,d} out of {counters['total']:,d} loci in the Illumina 174k allele frequency catalog")
print(f"Annotated {counters['found_t2t_assemblies_histogram']:,d} out of {counters['total']:,d} loci in the T2T assemblies allele frequency catalog")
print(f"Wrote {counters['total']:,d} records to {args.output_path}")
#%%



# gs://str-truth-set-v2/filter_vcf/all_repeats_including_homopolymers/combined/joined_only_pure_repeats.78_samples.variants.tsv.gz