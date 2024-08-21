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
parser.add_argument("--force", action="store_true", help="Run annotation step even if the output files already exist")

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
	("known_disease_associated_loci", "https://raw.githubusercontent.com/broadinstitute/str-analysis/main/str_analysis/variant_catalogs/variant_catalog_without_offtargets.GRCh38.json"),
	("comprehensive_catalog_from_Chiu_et_al", "https://zenodo.org/records/11522276/files/hg38.v1.bed.gz"),
	("PopSTR_catalog", "https://storage.googleapis.com/str-truth-set/hg38/ref/other/popstr_catalog_v2.bed.gz"),
	("illumina_catalog", "https://storage.googleapis.com/str-truth-set/hg38/ref/other/illumina_variant_catalog.sorted.bed.gz"),
	("GangSTR_v17", "https://s3.amazonaws.com/gangstr/hg38/genomewide/hg38_ver17.bed.gz"),
	("HipSTR_catalog", "https://github.com/HipSTR-Tool/HipSTR-references/raw/master/human/hg38.hipstr_reference.bed.gz"),
	("perfect_repeats_in_hg38", "https://storage.googleapis.com/str-truth-set/hg38/ref/other/colab-repeat-finder/hg38_repeats.motifs_1_to_1000bp.repeats_3x_and_spans_9bp/hg38_repeats.motifs_1_to_1000bp.repeats_3x_and_spans_9bp.bed.gz"),
	("polymorphic_loci_in_HPRC_assemblies", "https://storage.googleapis.com/str-truth-set-v2/filter_vcf/all_repeats_including_homopolymers_keeping_loci_that_have_overlapping_variants/combined/combined.51_samples.positive_loci.json"),
	("platinumTRs_v1.0", "https://zenodo.org/records/13178746/files/human_GRCh38_no_alt_analysis_set.platinumTRs-v1.0.trgt.bed.gz"),
	("adotto_TRs_catalog_v1.2", "https://storage.googleapis.com/str-truth-set/hg38/ref/other/adotto_tr_catalog_v1.2.bed.gz"),
	("mukamel_VNTR_catalog", "https://storage.googleapis.com/str-truth-set/hg38/ref/other/mukamel_VNTR_catalog.bed.gz"),
	("vamos_catalog_v2.1", "https://storage.googleapis.com/str-truth-set/hg38/ref/other/vamos_catalog.v2.1.bed.gz"),
]


catalog_paths = {}
for catalog_name, url in catalogs_in_order:
	if not os.path.isfile(os.path.basename(url)):
		run(f"wget -O {os.path.basename(url)}.tmp -qnc {url} && mv {os.path.basename(url)}.tmp {os.path.basename(url)}")
	catalog_paths[catalog_name] = os.path.abspath(os.path.basename(url))

# preprocess catalog of known disease-associated loci: split compound definitions
output_path = catalog_paths['known_disease_associated_loci'].replace(".json", ".split.json")
if not os.path.isfile(output_path):
	run(f"python3 -m str_analysis.split_adjacent_loci_in_expansion_hunter_catalog {catalog_paths['known_disease_associated_loci']}")
catalog_paths['known_disease_associated_loci'] = output_path

# change motif definition for the RFC1 locus from AARRG => AAAAG since our catalog doesn't currently support IUPAC codes
run(f"sed -i 's/AARRG/AAAAG/g' {catalog_paths['known_disease_associated_loci']}")

# compute stats for primary disease-associated loci
with open(catalog_paths["known_disease_associated_loci"]) as f:
	known_disease_associated_loci = json.load(f)

primary_disease_associated_loci = [
	x for x in known_disease_associated_loci if x["Diseases"] and (
		x["LocusId"].startswith("HOXA") or x["LocusId"].startswith("ARX") or "_" not in x["LocusId"]
	)
]
assert len(primary_disease_associated_loci) == 63

primary_disease_associated_loci_path = catalog_paths["known_disease_associated_loci"].replace(
	".json", ".primary_disease_associated_loci.json")

with open(primary_disease_associated_loci_path, "wt") as f:
	json.dump(primary_disease_associated_loci, f, indent=4)

catalogs_in_order = [("primary_known_disease_associated_loci", primary_disease_associated_loci_path)] + catalogs_in_order
catalog_paths["primary_known_disease_associated_loci"] = primary_disease_associated_loci_path

# convert PlatinumTR catalog to ExpansionHunter catalog format for comparison
if "platinumTRs_v1.0" in catalog_paths:
	path_after_conversion = catalog_paths["platinumTRs_v1.0"].replace(".bed.gz", ".catalog.json.gz")
	if not os.path.isfile(path_after_conversion):
		run(f"python3 -u -m str_analysis.convert_trgt_catalog_to_expansion_hunter_catalog -r {args.hg38_reference_fasta} {catalog_paths['platinumTRs_v1.0']} -o {path_after_conversion}")
	catalog_paths["platinumTRs_v1.0"] = path_after_conversion
else:
	print("WARNING: PlatinumTRs catalog not included in list of catalogs")

# convert GangSTR catalog to ExpansionHunter catalog format for comparison
if "GangSTR_v17" in catalog_paths:
	path_after_conversion = catalog_paths["GangSTR_v17"].replace(".bed.gz", ".catalog.json.gz")
	if not os.path.isfile(path_after_conversion):
		run(f"python3 -u -m str_analysis.convert_gangstr_spec_to_expansion_hunter_variant_catalog --verbose {catalog_paths['GangSTR_v17']} -o {path_after_conversion}")
	catalog_paths["GangSTR_v17"] = path_after_conversion
else:
	print("WARNING: GangSTR catalog not included in list of catalogs")

# convert HipSTR catalog to ExpansionHunter catalog format for comparison
if "HipSTR_catalog" in catalog_paths:
	path_after_conversion = catalog_paths["HipSTR_catalog"].replace(".bed.gz", ".catalog.bed.gz")
	if not os.path.isfile(path_after_conversion):
		run(f"python3 -u {scripts_dir}/convert_hipstr_catalog_to_regular_bed_file.py {catalog_paths['HipSTR_catalog']} -o {path_after_conversion}")
	catalog_paths["HipSTR_catalog"] = f"{path_after_conversion}"
else:
	print("WARNING: HipSTR catalog not included in list of catalogs")

# convert Chiu et al catalog to ExpansionHunter catalog format for comparison
if "comprehensive_catalog_from_Chiu_et_al" in catalog_paths:
	output_path = "hg38.Chiu_et_al.v1.bed.gz"
	if not os.path.isfile(output_path):
		# convert the comprehensive catalog from Chiu et al. to a regular bed file (removing the header and converting to 0-basd coords)
		run(f"gunzip -c {catalog_paths['comprehensive_catalog_from_Chiu_et_al']} | tail -n +3 | cut -f 1-4 | awk 'BEGIN {{OFS=\"\\t\"}} {{ print( $1, $2 - 1, $3, $4 ) }}' | gzip -c - > {output_path}")
	catalog_paths["comprehensive_catalog_from_Chiu_et_al"] = output_path

all_stats_tsv_paths = {}
for catalog_name, path in catalog_paths.items():
	annotated_catalog_path = re.sub("(.json|.bed)(.gz)?$", "", path) + ".annotated.json.gz"
	if not os.path.isfile(annotated_catalog_path) or args.force:
		run(f"""python3 -m str_analysis.annotate_and_filter_str_catalog \
			--reference-fasta {args.hg38_reference_fasta} \
			--skip-gene-annotations \
			--skip-disease-loci-annotations \
			--output-path {annotated_catalog_path} \
			{path}""")

	stats_tsv_path = re.sub("(.json|.bed)(.gz)?$", "", annotated_catalog_path) + ".catalog_stats.tsv"
	if True or not os.path.isfile(stats_tsv_path):
		print(f"Generating {stats_tsv_path}")
		run(f"python3 -m str_analysis.compute_catalog_stats --verbose {annotated_catalog_path}")
	all_stats_tsv_paths[catalog_name] = stats_tsv_path

print(f"Combining stats from all {len(all_stats_tsv_paths)} catalogs")
dfs = []
for catalog_name, stats_tsv_path in all_stats_tsv_paths.items():
	df = pd.read_table(stats_tsv_path)
	df.rename(columns={"catalog": "filename"}, inplace=True)
	df["catalog"] = catalog_name
	dfs.append(df)
df = pd.concat(dfs)
print(f"Sorting {len(df)} rows")
df["total"] = df["total"].apply(lambda s: str(s).replace(",", "")).astype(int)
df.sort_values("total", inplace=True)
df.to_csv(f"combined_catalog_stats.all_{len(all_stats_tsv_paths)}_catalogs.tsv", sep="\t", index=False)


diff = time.time() - start_time
print(f"Done with comparisons. Took {diff//3600:.0f}h, {(diff%3600)//60:.0f}m, {diff%60:.0f}s")

