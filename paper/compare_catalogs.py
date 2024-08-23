import argparse
import datetime
import json
import os
import pandas as pd
import re
import subprocess
import time


parser = argparse.ArgumentParser()
parser.add_argument("--hg38-reference-fasta", default="hg38.fa", help="Path of hg38 reference genome FASTA file")
parser.add_argument("--dry-run", action="store_true", help="Print commands without running them")
parser.add_argument("--outer-join", action="store_true", help="Perform an outer join on the catalogs")
parser.add_argument("--force", action="store_true", help="Run annotation step even if the output files already exist")
parser.add_argument("-k", "--keyword", help="Only process catalogs that contain this keyword")
args = parser.parse_args()


def run(command):
	print("-"*100)
	command = re.sub("[ \\t]{2,}", "  ", command)  # remove extra spaces
	print(command)
	if not args.dry_run or command.startswith("mkdir"):
		subprocess.run(command, shell=True, check=True)

def chdir(d):
	print(f"cd {d}")
	os.chdir(d)

for key in "hg38_reference_fasta", :
	path = getattr(args, key)
	if not os.path.isfile(path):
		parser.error(f"{key} file not found {path}")
	setattr(args, key, os.path.abspath(path))

base_dir = os.path.abspath(".")
scripts_dir = os.path.abspath(f"../scripts")
working_dir = os.path.abspath(f"compare_catalogs")


run(f"mkdir -p {working_dir}")
chdir(working_dir)


start_time = time.time()

#  Download [f"https://github.com/DecodeGenetics/popSTR/raw/master/chr{i}markerInfo.gz" for i in range(1, 23)]

# list of source catalogs, in order. If a locus is defined in more than one catalog (ie. overlapping boundaries,
# same motif), then the definition in the catalog that's earlier in the list will take precedence over definitions in
# subsequent catalogs.
catalogs_in_order = [
	("KnownDiseaseAssociatedLoci", "https://raw.githubusercontent.com/broadinstitute/str-analysis/main/str_analysis/variant_catalogs/variant_catalog_without_offtargets.GRCh38.json"),
	("Illumina174kPolymorphicTRs", "https://storage.googleapis.com/str-truth-set/hg38/ref/other/illumina_variant_catalog.sorted.bed.gz"),
	("UCSC_SimpleRepeatTrack", "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/simpleRepeat.txt.gz"),
	("VamosCatalog_v2.1", "https://storage.googleapis.com/str-truth-set/hg38/ref/other/vamos_catalog.v2.1.bed.gz"),
	("GangSTR_v17", "https://s3.amazonaws.com/gangstr/hg38/genomewide/hg38_ver17.bed.gz"),
	("HipSTR_Catalog", "https://github.com/HipSTR-Tool/HipSTR-references/raw/master/human/hg38.hipstr_reference.bed.gz"),
	("Adotto_v1.2", "https://storage.googleapis.com/str-truth-set/hg38/ref/other/adotto_tr_catalog_v1.2.bed.gz"),
	("PopSTR_Catalog", "https://storage.googleapis.com/str-truth-set/hg38/ref/other/popstr_catalog_v2.bed.gz"),
	("PlatinumTRs_v1.0", "https://zenodo.org/records/13178746/files/human_GRCh38_no_alt_analysis_set.platinumTRs-v1.0.trgt.bed.gz"),
	("Chiu_et_al", "https://zenodo.org/records/11522276/files/hg38.v1.bed.gz"),
	("PerfectRepeatsInReference", "https://storage.googleapis.com/str-truth-set/hg38/ref/other/colab-repeat-finder/hg38_repeats.motifs_1_to_1000bp.repeats_3x_and_spans_9bp/hg38_repeats.motifs_1_to_1000bp.repeats_3x_and_spans_9bp.bed.gz"),
	("PolymorphicTRsInT2TAssemblies", "https://storage.googleapis.com/str-truth-set-v2/filter_vcf/all_repeats_including_homopolymers_keeping_loci_that_have_overlapping_variants/combined/combined.78_samples.positive_loci.json"),
	("MukamelVNTRs", "https://storage.googleapis.com/str-truth-set/hg38/ref/other/mukamel_VNTR_catalog.bed.gz"),
]


catalog_paths = {}
for catalog_name, url in catalogs_in_order:
	if not os.path.isfile(os.path.basename(url)):
		run(f"wget -O {os.path.basename(url)}.tmp -qnc {url} && mv {os.path.basename(url)}.tmp {os.path.basename(url)}")
	catalog_paths[catalog_name] = os.path.abspath(os.path.basename(url))

# preprocess catalog of known disease-associated loci: split compound definitions
output_path = catalog_paths['KnownDiseaseAssociatedLoci'].replace(".json", ".split.json")
if not os.path.isfile(output_path):
	run(f"python3 -m str_analysis.split_adjacent_loci_in_expansion_hunter_catalog {catalog_paths['KnownDiseaseAssociatedLoci']}")
catalog_paths['KnownDiseaseAssociatedLoci'] = output_path

# change motif definition for the RFC1 locus from AARRG => AAAAG since our catalog doesn't currently support IUPAC codes
run(f"sed -i 's/AARRG/AAAAG/g' {catalog_paths['KnownDiseaseAssociatedLoci']}")

# compute stats for primary disease-associated loci
with open(catalog_paths["KnownDiseaseAssociatedLoci"]) as f:
	known_disease_associated_loci = json.load(f)

primary_disease_associated_loci = [
	x for x in known_disease_associated_loci if x["Diseases"] and (
		x["LocusId"].startswith("HOXA") or x["LocusId"].startswith("ARX") or "_" not in x["LocusId"]
	)
]
assert len(primary_disease_associated_loci) == 63

primary_disease_associated_loci_path = catalog_paths["KnownDiseaseAssociatedLoci"].replace(
	".json", ".primary_disease_associated_loci.json")

with open(primary_disease_associated_loci_path, "wt") as f:
	json.dump(primary_disease_associated_loci, f, indent=4)

catalogs_in_order = [("PrimaryKnownDiseaseAssociatedLoci", primary_disease_associated_loci_path)] + catalogs_in_order
catalog_paths["PrimaryKnownDiseaseAssociatedLoci"] = primary_disease_associated_loci_path

# convert PlatinumTR catalog to ExpansionHunter catalog format for comparison
if "PlatinumTRs_v1.0" in catalog_paths and (
	not args.keyword or args.keyword.lower() in "platinumTRs_v1.0".lower() or args.keyword.lower() in catalog_paths["platinumTRs_v1.0"].lower()):
	path_after_conversion = catalog_paths["platinumTRs_v1.0"].replace(".bed.gz", ".catalog.json.gz")
	if not os.path.isfile(path_after_conversion):
		run(f"python3 -u -m str_analysis.convert_trgt_catalog_to_expansion_hunter_catalog -r {args.hg38_reference_fasta} {catalog_paths['platinumTRs_v1.0']} -o {path_after_conversion}")
	catalog_paths["platinumTRs_v1.0"] = path_after_conversion
else:
	print("WARNING: PlatinumTRs catalog not included in list of catalogs")

# convert GangSTR catalog to ExpansionHunter catalog format for comparison
if "GangSTR_v17" in catalog_paths and (
	not args.keyword or args.keyword.lower() in "GangSTR_v17".lower() or args.keyword.lower() in catalog_paths["GangSTR_v17"].lower()):
	path_after_conversion = catalog_paths["GangSTR_v17"].replace(".bed.gz", ".catalog.json.gz")
	if not os.path.isfile(path_after_conversion):
		run(f"python3 -u -m str_analysis.convert_gangstr_spec_to_expansion_hunter_variant_catalog --verbose {catalog_paths['GangSTR_v17']} -o {path_after_conversion}")
	catalog_paths["GangSTR_v17"] = path_after_conversion
else:
	print("WARNING: GangSTR catalog not included in list of catalogs")

# convert HipSTR catalog to ExpansionHunter catalog format for comparison
if "HipSTR_catalog" in catalog_paths and (
	not args.keyword or args.keyword.lower() in "HipSTR_catalog".lower() or args.keyword.lower() in catalog_paths["HipSTR_catalog"].lower()):
	path_after_conversion = catalog_paths["HipSTR_catalog"].replace(".bed.gz", ".catalog.bed.gz")
	if not os.path.isfile(path_after_conversion):
		run(f"python3 -u {scripts_dir}/convert_hipstr_catalog_to_regular_bed_file.py {catalog_paths['HipSTR_catalog']} -o {path_after_conversion}")
	catalog_paths["HipSTR_catalog"] = path_after_conversion
else:
	print("WARNING: HipSTR catalog not included in list of catalogs")

# convert UCSC track to BED format
if "UCSC_SimpleRepeatTrack" in catalog_paths and (
	not args.keyword or args.keyword.lower() in "UCSC_SimpleRepeatTrack".lower() or args.keyword.lower() in catalog_paths["UCSC_SimpleRepeatTrack"].lower()):
	path_after_conversion = catalog_paths["UCSC_SimpleRepeatTrack"].replace(".txt.gz", "") + "_track_from_UCSC.bed.gz"
	if not os.path.isfile(path_after_conversion):
		run(f"python3 -u {scripts_dir}/convert_UCSC_SimpleRepeatTrack_to_bed.py {catalog_paths['UCSC_SimpleRepeatTrack']} -o {path_after_conversion}")
	catalog_paths["UCSC_SimpleRepeatTrack"] = path_after_conversion

# convert Chiu et al catalog to ExpansionHunter catalog format for comparison
if "Chiu_et_al" in catalog_paths and (
	not args.keyword or args.keyword.lower() in "Chiu_et_al".lower() or args.keyword.lower() in catalog_paths["Chiu_et_al"].lower()):
	output_path = "hg38.Chiu_et_al.v1.bed.gz"
	if not os.path.isfile(output_path):
		# convert the comprehensive catalog from Chiu et al. to a regular bed file (removing the header and converting to 0-basd coords)
		run(f"gunzip -c {catalog_paths['Chiu_et_al']} | tail -n +3 | cut -f 1-4 | awk 'BEGIN {{OFS=\"\\t\"}} {{ print( $1, $2 - 1, $3, $4 ) }}' | gzip -c - > {output_path}")
	catalog_paths["Chiu_et_al"] = output_path

all_stats_tsv_paths = {}
catalog_paths_for_outer_join = []
for catalog_name, _ in catalogs_in_order:
	path = catalog_paths[catalog_name]

	if args.keyword and args.keyword.lower() not in catalog_name.lower() and args.keyword.lower() not in path.lower():
		continue

	annotated_catalog_path = re.sub("(.json|.bed)(.gz)?$", "", path) + ".annotated.json.gz"
	if not os.path.isfile(annotated_catalog_path) or args.force:
		run(f"""python3 -m str_analysis.annotate_and_filter_str_catalog \
			--verbose \
			--trim-loci \
			--reference-fasta {args.hg38_reference_fasta} \
			--skip-gene-annotations \
			--skip-disease-loci-annotations \
			--output-path {annotated_catalog_path} \
			{path}""")

	if args.outer_join:
		catalog_paths_for_outer_join.append(f"{catalog_name.title().replace(' ', '_')}:{annotated_catalog_path}")

	# just compute stats
	stats_tsv_path = re.sub("(.json|.bed)(.gz)?$", "", annotated_catalog_path) + ".catalog_stats.tsv"
	if not os.path.isfile(stats_tsv_path) or args.force:
		print(f"Generating {stats_tsv_path}")
		run(f"python3 -m str_analysis.compute_catalog_stats --verbose {annotated_catalog_path}")

	all_stats_tsv_paths[catalog_name] = stats_tsv_path



if args.outer_join:
	run(f"""python3 -u -m str_analysis.merge_loci \
		--add-source-field \
		--add-found-in-fields \
		--output-format JSON \
		--overlapping-loci-action keep-first \
		--write-merge-stats-tsv \
		--write-outer-join-overlap-table \
		--outer-join-overlap-table-min-sources 2 \
		--output-prefix merged.all_{len(all_stats_tsv_paths)}_catalogs \
		{' '.join(catalog_paths_for_outer_join)}""")

print(f"Combining stats from all {len(all_stats_tsv_paths)} catalogs")
dfs = []
for catalog_name, stats_tsv_path in all_stats_tsv_paths.items():
	df = pd.read_table(stats_tsv_path)
	df.rename(columns={"catalog": "filename"}, inplace=True)
	df["catalog"] = catalog_name
	dfs.append(df)

if dfs:
	df = pd.concat(dfs)
	print(f"Sorting {len(df)} rows")
	df["total"] = df["total"].apply(lambda s: str(s).replace(",", "")).astype(int)
	df.sort_values("total", inplace=True)
	df.to_csv(f"combined_catalog_stats.all_{len(all_stats_tsv_paths)}_catalogs.tsv", sep="\t", index=False)
else:
	print("No catalogs to combine")

diff = time.time() - start_time
print(f"Done. Took {diff//3600:.0f}h, {(diff%3600)//60:.0f}m, {diff%60:.0f}s")

