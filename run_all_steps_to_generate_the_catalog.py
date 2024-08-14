import argparse
import datetime
import json
import os
import re
import subprocess
import time

# install str-analysis python package
os.system("""
if ! python3 -m pip show str-analysis &> /dev/null
then
    python3 -m pip install --upgrade git+https://github.com/broadinstitute/str-analysis.git
fi
""")


parser = argparse.ArgumentParser()
parser.add_argument("--hg38-reference-fasta", default="hg38.fa", help="Path of hg38 reference genome FASTA file")
parser.add_argument("--gencode-gtf", default="gencode.v46.basic.annotation.gtf.gz", help="Gene annotations GTF file")
parser.add_argument("--variation-clusters-bed", default="vcs_merged.bed.gz", help="Variation clusters file shared by Egor Dolzhenko")
parser.add_argument("--skip-variation-cluster-annotations", action="store_true", help="Skip adding variation cluster annotations to the catalog")
parser.add_argument("--output-prefix", default="repeat_catalog_v1.hg38")
parser.add_argument("--dry-run", action="store_true", help="Print commands without running them")

args = parser.parse_args()


def run(command):
	command = re.sub("[ \\t]{2,}", "  ", command)  # remove extra spaces
	print(command)
	if not args.dry_run or command.startswith("mkdir"):
		subprocess.run(command, shell=True, check=True)

def chdir(d):
	print(f"cd {d}")
	os.chdir(d)

for key in "hg38_reference_fasta", "gencode_gtf", "variation_clusters_bed":
	path = getattr(args, key)
	if not os.path.isfile(path):
		if key == "variation_clusters_bed" and args.skip_variation_cluster_annotations:
			setattr(args, key, None)
		else:
			parser.error(f"{key} file not found {path}")
	setattr(args, key, os.path.abspath(path))

base_dir = os.path.abspath(".")
working_dir = os.path.abspath(f"results__{datetime.datetime.now().strftime('%F')}")

run(f"mkdir -p {working_dir}")
chdir(working_dir)

# create a release draft folder
timestamp = datetime.datetime.now().strftime('%Y-%m-%d')
release_draft_folder = os.path.abspath(f"release_draft_{timestamp}")
run(f"mkdir -p {release_draft_folder}")


start_time = time.time()

# STEP #1  (already done for hg38, and results are publicly available. Will just download them in step #2)
#git clone git@github.com:broadinstitute/colab-repeat-finder.git
#cd colab-repeat-finder/python
#python3 perfect_repeat_finder.py  --min-repeats 3  --min-span 9  --min-motif-size 1  --max-motif-size 1000  --output-prefix perfect_repeats.hg38  --show-progress-bar  ${REFERENCE_FASTA_PATH}


# STEP #2:  generate catalog

# list of source catalogs, in order. If a locus is defined in more than one catalog (ie. overlapping boundaries,
# same motif), then the definition in the catalog that's earlier in the list will take precedence over definitions in
# subsequent catalogs.
source_catalogs_in_order = [
	("KNOWN_DISEASE_ASSOCIATED_LOCI", "https://raw.githubusercontent.com/broadinstitute/str-analysis/main/str_analysis/variant_catalogs/variant_catalog_without_offtargets.GRCh38.json"),
	("ILLUMINA_CATALOG", "https://storage.googleapis.com/str-truth-set/hg38/ref/other/illumina_variant_catalog.sorted.bed.gz"),
	("ALL_PERFECT_REPEATS_CATALOG", "https://storage.googleapis.com/str-truth-set/hg38/ref/other/colab-repeat-finder/hg38_repeats.motifs_1_to_1000bp.repeats_3x_and_spans_9bp/hg38_repeats.motifs_1_to_1000bp.repeats_3x_and_spans_9bp.bed.gz"),
	("TRUTH_SET_CATALOG", "https://storage.googleapis.com/str-truth-set-v2/filter_vcf/all_repeats_including_homopolymers_keeping_loci_that_have_overlapping_variants/combined/combined.51_samples.positive_loci.json")
]


source_catalog_paths = {}
for catalog_name, url in source_catalogs_in_order:
	run(f"wget -qnc {url}")
	source_catalog_paths[catalog_name] = os.path.abspath(os.path.basename(url))

# preprocess catalog of known disease-associated loci: split compound definitions
run(f"python3 -m str_analysis.split_adjacent_loci_in_expansion_hunter_catalog {source_catalog_paths['KNOWN_DISEASE_ASSOCIATED_LOCI']}")
source_catalog_paths['KNOWN_DISEASE_ASSOCIATED_LOCI'] = source_catalog_paths['KNOWN_DISEASE_ASSOCIATED_LOCI'].replace(".json", ".split.json")
# change motif definition for the RFC1 locus from AARRG => AAAAG since our catalog doesn't currently support IUPAC codes
run(f"sed -i 's/AARRG/AAAAG/g' {source_catalog_paths['KNOWN_DISEASE_ASSOCIATED_LOCI']}")


# compute stats for primary disease-associated loci
with open(source_catalog_paths["KNOWN_DISEASE_ASSOCIATED_LOCI"]) as f:
	known_disease_associated_loci = json.load(f)

primary_disease_associated_loci = [
	x for x in known_disease_associated_loci if x["Diseases"] and (
		x["LocusId"].startswith("HOXA") or x["LocusId"].startswith("ARX") or "_" not in x["LocusId"]
	)
]
# check the number of primary disease-associated loci, not counting adjacent repeats and historic candidate loci that
# are not currently considered monogenic
assert len(primary_disease_associated_loci) == 63

primary_disease_associated_loci_path = source_catalog_paths["KNOWN_DISEASE_ASSOCIATED_LOCI"].replace(
	".json", ".primary_disease_associated_loci.json")

with open(primary_disease_associated_loci_path, "wt") as f:
	json.dump(primary_disease_associated_loci, f, indent=4)

run(f"""python3 -u -m str_analysis.annotate_and_filter_str_catalog \
	--verbose \
	--reference-fasta {args.hg38_reference_fasta} \
	--min-interval-size-bp 1 \
	--skip-gene-annotations \
	--skip-mappability-annotations \
	--skip-disease-loci-annotations \
	--discard-loci-with-non-ACGT-bases-in-reference \
	--discard-loci-with-non-ACGTN-bases-in-motif \
	--output-path {primary_disease_associated_loci_path} \
	{primary_disease_associated_loci_path}""")

run(f"python3 -m str_analysis.compute_catalog_stats --verbose {primary_disease_associated_loci_path}")

for motif_size_label, min_motif_size, max_motif_size in [
	("1_to_1000bp",  1, 1000),
	("2_to_1000bp",  2, 1000),
	("homopolymers", 1, 1),
	("2_to_6bp",     2, 6),
	("7_to_1000bp",  7, 1000)
]:
	print("="*200)
	chdir(working_dir)
	run(f"mkdir -p {motif_size_label}")
	chdir(motif_size_label)

	output_prefix = os.path.abspath(f"{args.output_prefix}.{motif_size_label}")

	filtered_source_catalog_paths = {}
	for catalog_name, _ in source_catalogs_in_order:
		catalog_path = source_catalog_paths[catalog_name]
		if catalog_path.endswith(".json") or catalog_path.endswith(".json.gz"):
			filtered_catalog_path = catalog_path.replace(".json", ".filtered.json.gz")
		elif catalog_path.endswith(".bed.gz"):
			filtered_catalog_path = catalog_path.replace(".bed.gz", ".filtered.json.gz")
		else:
			raise ValueError(f"Unexpected file extension for {catalog_path}")

		filtered_catalog_path = os.path.abspath(os.path.basename(filtered_catalog_path))
		filtered_source_catalog_paths[catalog_name] = filtered_catalog_path

		run(f"""python3 -u -m str_analysis.annotate_and_filter_str_catalog \
			--verbose \
			--reference-fasta {args.hg38_reference_fasta} \
			--min-motif-size {min_motif_size} \
			--max-motif-size {max_motif_size} \
			--min-interval-size-bp 1 \
			--skip-gene-annotations \
			--skip-mappability-annotations \
			--skip-disease-loci-annotations \
			--discard-loci-with-non-ACGT-bases-in-reference \
			--discard-loci-with-non-ACGTN-bases-in-motif \
			--output-path {filtered_catalog_path} \
			{catalog_path}""")

		print(f"Stats for {catalog_path}")
		run(f"python3 -m str_analysis.compute_catalog_stats --verbose {filtered_catalog_path}")

	catalog_paths = " ".join([filtered_source_catalog_paths[catalog_name] for catalog_name, _ in source_catalogs_in_order])
	run(f"""python3 -u -m str_analysis.merge_loci \
		--merge-adjacent-loci-with-same-motif \
		--add-source-field \
		--output-format JSON \
		--write-bed-files-with-new-loci \
		--output-prefix {output_prefix}.merged \
		{catalog_paths}""")

	run(f"""python3 -u -m str_analysis.annotate_and_filter_str_catalog \
		--reference-fasta {args.hg38_reference_fasta} \
		--genes-gtf {args.gencode_gtf} \
		--gene-models-source gencode \
		--min-motif-size {min_motif_size} \
		--max-motif-size {max_motif_size} \
		--min-interval-size-bp 1 \
		--discard-overlapping-intervals-with-similar-motifs \
		--output-path {output_prefix}.merged_and_annotated.json.gz \
		{output_prefix}.merged.json.gz""")


	# Replace "Source" field filenames with more user-friendly source names
	run(f"""cp {output_prefix}.merged_and_annotated.json.gz  temp.json.gz
	gunzip -c temp.json.gz | sed 's/{os.path.basename(filtered_source_catalog_paths['KNOWN_DISEASE_ASSOCIATED_LOCI'])}/known disease-associated loci/' | gzip -c - > temp2.json.gz && mv temp2.json.gz temp.json.gz
	gunzip -c temp.json.gz | sed 's/{os.path.basename(filtered_source_catalog_paths['ILLUMINA_CATALOG'])}/Illumina catalog of 174k polymorphic loci/'  | gzip -c - > temp2.json.gz && mv temp2.json.gz temp.json.gz
	gunzip -c temp.json.gz | sed 's/{os.path.basename(filtered_source_catalog_paths['ALL_PERFECT_REPEATS_CATALOG'])}/perfect repeats in hg38/'  | gzip -c - > temp2.json.gz && mv temp2.json.gz temp.json.gz
	gunzip -c temp.json.gz | sed 's/{os.path.basename(filtered_source_catalog_paths['TRUTH_SET_CATALOG'])}/polymorphic TR loci in HPRC assemblies/'  | gzip -c - > temp2.json.gz && mv temp2.json.gz temp.json.gz
	mv temp.json.gz {output_prefix}.merged_and_annotated.json.gz""")

	if args.variation_clusters_bed and motif_size_label != "homopolymers":
		run(f"""python3 {base_dir}/scripts/add_variation_cluster_annotations_to_catalog.py \
			--verbose \
			--output-json-path {output_prefix}.merged_and_annotated.with_variation_clusters.json.gz \
			{args.variation_clusters_bed} \
			{output_prefix}.merged_and_annotated.json.gz""")

		run(f"mv {output_prefix}.merged_and_annotated.with_variation_clusters.json.gz {output_prefix}.merged_and_annotated.json.gz")

	# STEP #3: convert the catalog from ExpansionHunter catalog format to bed, TRGT, LongTR, HipSTR, and GangSTR formats
	run(f"python3 -m str_analysis.convert_expansion_hunter_variant_catalog_to_bed --split-adjacent-repeats {output_prefix}.merged_and_annotated.json.gz  --output-file {output_prefix}.bed.gz")
	run(f"python3 -m str_analysis.convert_expansion_hunter_variant_catalog_to_trgt_catalog --split-adjacent-repeats {output_prefix}.merged_and_annotated.json.gz  --output-file {output_prefix}.TRGT.bed")
	run(f"python3 -m str_analysis.convert_expansion_hunter_variant_catalog_to_longtr_format  {output_prefix}.merged_and_annotated.json.gz  --output-file {output_prefix}.LongTR.bed")
	run(f"python3 -m str_analysis.convert_expansion_hunter_variant_catalog_to_hipstr_format  {output_prefix}.merged_and_annotated.json.gz  --output-file {output_prefix}.HipSTR.bed")
	run(f"python3 -m str_analysis.convert_expansion_hunter_variant_catalog_to_gangstr_spec   {output_prefix}.merged_and_annotated.json.gz  --output-file {output_prefix}.GangSTR.bed")

	# Confirm that the TRGT catalog passes 'trgt validate'
	run(f"trgt validate --genome {args.hg38_reference_fasta}  --repeats {output_prefix}.TRGT.bed")

	# Spot-check several known disease-associated loci to make sure they made it into the final catalog
	#if min_motif_size <= 3 and max_motif_size >= 24:
	#	for end_coord in 3074933, 69037304, 39348479, 41746032, 80147139:    # HTT, FXN, RFC1, PHOX2B, EIF4A3
	#		run(f"grep $'\t{end_coord}$'  {output_prefix}.TRGT.bed")

	# STEP #4: copy files to the release_draft folder and compute catalog stats
	for path in (f"{output_prefix}.bed.gz",
				 f"{output_prefix}.bed.gz.tbi",
				 f"{output_prefix}.merged_and_annotated.json.gz",
				 f"{output_prefix}.TRGT.bed",
				 f"{output_prefix}.LongTR.bed",
				 f"{output_prefix}.HipSTR.bed",
				 f"{output_prefix}.GangSTR.bed"):
		run(f"cp {path} {release_draft_folder}")

	run(f"python3 -m str_analysis.compute_catalog_stats --verbose {output_prefix}.merged_and_annotated.json.gz")

	# report hours, minutes, seconds relative to start_time
	diff = time.time() - start_time
	print(f"Done generating {output_prefix} catalog. Took {diff//3600:.0f}h, {(diff%3600)//60:.0f}m, {diff%60:.0f}s")

	if motif_size_label != "1_to_1000bp":
		continue

	# compare to the GangSTR_v17 catalog to make sure it's included
	comparison_catalogs_in_order = [
		("GangSTR_v17", "https://s3.amazonaws.com/gangstr/hg38/genomewide/hg38_ver17.bed.gz"),
	]

	comparison_catalog_paths = {}
	for catalog_name, url in comparison_catalogs_in_order:
		run(f"wget -qnc {url}")
		comparison_catalog_paths[catalog_name] = os.path.abspath(os.path.basename(url))

	path_after_conversion = comparison_catalog_paths["GangSTR_v17"].replace(".bed.gz", ".json.gz")
	run(f"python3 -u -m str_analysis.convert_gangstr_spec_to_expansion_hunter_variant_catalog --verbose {comparison_catalog_paths['GangSTR_v17']} -o {path_after_conversion}")
	comparison_catalog_paths["GangSTR_v17"] = path_after_conversion

	# STEP #5:  compare catalog to other catalogs
	start_time = time.time()

	for catalog_name, path in comparison_catalog_paths.items():
		filtered_comparison_catalog_path = re.sub("(.json.gz|.bed.gz)$", ".filtered.json.gz", path)

		run(f"""python3 -m str_analysis.annotate_and_filter_str_catalog \
			--reference-fasta {args.hg38_reference_fasta} \
			--skip-gene-annotations \
			--skip-mappability-annotations \
			--skip-disease-loci-annotations \
			--min-motif-size {min_motif_size} \
			--max-motif-size {max_motif_size} \
			--output-path {filtered_comparison_catalog_path} \
			--verbose \
			{path}""")

		run(f"python3 -m str_analysis.compute_catalog_stats --verbose {filtered_comparison_catalog_path}")

		run(f"""python3 -u -m str_analysis.merge_loci \
			--verbose \
			--output-prefix {catalog_name} \
			--write-bed-files-with-new-loci \
			{output_prefix}.merged_and_annotated.json.gz \
			{filtered_comparison_catalog_path}""")

	diff = time.time() - start_time
	print(f"Done with comparisons. Took {diff//3600:.0f}h, {(diff%3600)//60:.0f}m, {diff%60:.0f}s")

