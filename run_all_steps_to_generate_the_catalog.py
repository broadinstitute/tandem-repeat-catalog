import argparse
import datetime
import gzip
import json
import os
import re
import subprocess
import time

from str_analysis.utils.file_utils import file_exists

# install str-analysis python package
os.system("""
if ! python3 -m pip show str-analysis &> /dev/null
then
	python3 -m pip install --upgrade git+https://github.com/broadinstitute/str-analysis.git
fi
""")


def run(command, step_number=None):
	command = re.sub("[ \\t]{2,}", "  ", command)  # remove extra spaces
	if not args.dry_run or command.startswith("mkdir"):
		if step_number is None or (
		   	not (args.only_step is not None and step_number != args.only_step) and
			not (args.start_with_step is not None and step_number < args.start_with_step) and
			not (args.end_with_step is not None and step_number > args.end_with_step)
		):
			if step_number is not None:
				print(f"STEP #{step_number}: {command}")
			else:
				print(command)
			subprocess.run(command, shell=True, check=True)

def chdir(d):
	print(f"cd {d}")
	os.chdir(d)



parser = argparse.ArgumentParser()
parser.add_argument("--hg38-reference-fasta", default="hg38.fa", help="Path of hg38 reference genome FASTA file")
parser.add_argument("--gencode-gtf", default="gencode.v46.basic.annotation.gtf.gz", help="Gene annotations GTF file")
parser.add_argument("--output-prefix", default="repeat_catalog_v1.hg38")
parser.add_argument("--only-step", type=int, help="Only run this one step")
parser.add_argument("--start-with-step", type=int, help="Start with a specific step number")
parser.add_argument("--end-with-step", type=int, help="End with a specific step number")
parser.add_argument("--lps-annotations", default="gs://tandem-repeat-catalog/v1.0/HPRC_100_LongestPureSegmentQuantiles.txt.gz",
					help="Path of the LPS annotations table shared by Matt Danzi")
parser.add_argument("--skip-lps-annotations", action="store_true",
					help="Skip adding variation cluster annotations to the catalog")
parser.add_argument("--variation-clusters-bed", default="gs://tandem-repeat-catalog/v1.0/vcs_v1.0.bed.gz",
					help="Variation clusters file shared by Egor Dolzhenko")
parser.add_argument("--skip-variation-cluster-annotations", action="store_true",
					help="Skip adding variation cluster annotations to the catalog")
parser.add_argument("--variation-clusters-output-prefix", default="variation_clusters_v1.hg38")
parser.add_argument("--timestamp", default=datetime.datetime.now().strftime('%Y-%m-%d'),
					help="Timestamp to use in the output directory name")
parser.add_argument("--dry-run", action="store_true", help="Print commands without running them")

args = parser.parse_args()

print("TIMESTAMP:", args.timestamp)


for key in "hg38_reference_fasta", "gencode_gtf", "variation_clusters_bed", "lps_annotations":
	if (
		key == "lps_annotations" and args.skip_lps_annotations
	) or (
		key == "variation_clusters_bed" and args.skip_variation_cluster_annotations
	):
		setattr(args, key, None)
		continue
	path = getattr(args, key)
	if not file_exists(path):
		parser.error(f"{key} file not found {path}")

	if path.startswith("gs://"):
		run(f"gsutil -m cp {path} .")
		path = os.path.basename(path)

	setattr(args, key, os.path.abspath(path))

base_dir = os.path.abspath(".")
working_dir = os.path.abspath(f"results__{args.timestamp}")

run(f"mkdir -p {working_dir}")
chdir(working_dir)

# create a release draft folder
release_draft_folder = os.path.abspath(f"release_draft_{args.timestamp}")
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
	("KnownDiseaseAssociatedLoci", "https://raw.githubusercontent.com/broadinstitute/str-analysis/main/str_analysis/variant_catalogs/variant_catalog_without_offtargets.GRCh38.json"),
	("Illumina174kPolymorphicTRs", "https://storage.googleapis.com/str-truth-set/hg38/ref/other/illumina_variant_catalog.sorted.bed.gz"),
	("PerfectRepeatsInReference", "https://storage.googleapis.com/str-truth-set/hg38/ref/other/colab-repeat-finder/hg38_repeats.motifs_1_to_1000bp.repeats_3x_and_spans_9bp/hg38_repeats.motifs_1_to_1000bp.repeats_3x_and_spans_9bp.bed.gz"),
	("PolymorphicTRsInT2TAssemblies", "https://storage.googleapis.com/str-truth-set-v2/filter_vcf/all_repeats_including_homopolymers_keeping_loci_that_have_overlapping_variants/combined/merged_expansion_hunter_catalog.78_samples.json.gz")
]


source_catalog_paths = {}
for catalog_name, url in source_catalogs_in_order:
	if not os.path.isfile(os.path.basename(url)):
		run(f"wget -O {os.path.basename(url)}.tmp -qnc {url} && mv {os.path.basename(url)}.tmp {os.path.basename(url)}")
	source_catalog_paths[catalog_name] = os.path.abspath(os.path.basename(url))

# preprocess catalog of known disease-associated loci: split compound definitions
run(f"python3 -u -m str_analysis.split_adjacent_loci_in_expansion_hunter_catalog {source_catalog_paths['KnownDiseaseAssociatedLoci']}", step_number=0)
source_catalog_paths['KnownDiseaseAssociatedLoci'] = source_catalog_paths['KnownDiseaseAssociatedLoci'].replace(".json", ".split.json")
# change motif definition for the RFC1 locus from AARRG => AAAAG since our catalog doesn't currently support IUPAC codes
run(f"sed -i 's/AARRG/AAAAG/g' {source_catalog_paths['KnownDiseaseAssociatedLoci']}", step_number=0)
run(f"gzip -f {source_catalog_paths['KnownDiseaseAssociatedLoci']}", step_number=0)

source_catalog_paths['KnownDiseaseAssociatedLoci'] += ".gz"

# compute stats for primary disease-associated loci
if not args.dry_run:
	with gzip.open(source_catalog_paths["KnownDiseaseAssociatedLoci"]) as f:
		known_disease_associated_loci = json.load(f)

	primary_disease_associated_loci = [
		x for x in known_disease_associated_loci if x["Diseases"] and (
			x["LocusId"].startswith("HOXA") or x["LocusId"].startswith("ARX") or "_" not in x["LocusId"]
		)
	]
	# check the number of primary disease-associated loci, not counting adjacent repeats and historic candidate loci that
	# are not currently considered monogenic
	assert len(primary_disease_associated_loci) == 63

	primary_disease_associated_loci_path = source_catalog_paths["KnownDiseaseAssociatedLoci"].replace(
		".json.gz", ".primary_disease_associated_loci.json.gz")

	with gzip.open(primary_disease_associated_loci_path, "wt") as f:
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
		{primary_disease_associated_loci_path}""", step_number=1)

	run(f"python3 -m str_analysis.compute_catalog_stats --verbose {primary_disease_associated_loci_path}", step_number=2)

adjacent_repeats_source_bed = None
for motif_size_label, min_motif_size, max_motif_size, release_tar_gz_path in [
	("1_to_1000bp_motifs",  1, 1000, None),
	#("2_to_1000bp_motifs",  2, 1000, f"{args.output_prefix}.subset.all_loci_except_homopolymers.tar.gz"),
	#("homopolymers",        1, 1,    f"{args.output_prefix}.subset.only_homopolymer_loci.tar.gz"),
	#("2_to_6bp_motifs",     2, 6,    f"{args.output_prefix}.subset.only_loci_with_2_to_6bp_motifs.tar.gz"),
	#("7_to_1000bp_motifs",  7, 1000, f"{args.output_prefix}.subset.only_loci_with_7_to_1000bp_motifs.tar.gz"),
]:

	if (args.start_with_step and args.start_with_step > 2) and motif_size_label != "1_to_1000bp_motifs":
		continue


	print("="*200)
	chdir(working_dir)
	run(f"mkdir -p {motif_size_label}")
	chdir(motif_size_label)

	output_prefix = os.path.abspath(f"{args.output_prefix}.{motif_size_label}")

	filtered_source_catalog_paths = {}
	for catalog_name, _ in source_catalogs_in_order:
		catalog_path = source_catalog_paths[catalog_name]
		if catalog_path.endswith(".json.gz"):
			filtered_catalog_path = catalog_path.replace(".json.gz", ".filtered.json.gz")
		elif catalog_path.endswith(".bed.gz"):
			filtered_catalog_path = catalog_path.replace(".bed.gz", ".filtered.json.gz")
		else:
			raise ValueError(f"Unexpected file extension for {catalog_path}")

		filtered_catalog_path = os.path.abspath(os.path.basename(filtered_catalog_path))
		filtered_source_catalog_paths[catalog_name] = filtered_catalog_path

		run(f"""python3 -u -m str_analysis.annotate_and_filter_str_catalog \
			--verbose \
			--known-disease-associated-loci {primary_disease_associated_loci_path} \
			--reference-fasta {args.hg38_reference_fasta} \
			--min-motif-size {min_motif_size} \
			--max-motif-size {max_motif_size} \
			--min-interval-size-bp 1 \
			--skip-gene-annotations \
			--skip-mappability-annotations \
			--skip-disease-loci-annotations \
			--set-locus-id \
			--discard-loci-with-non-ACGT-bases-in-reference \
			--discard-loci-with-non-ACGTN-bases-in-motif \
			--output-path {filtered_catalog_path} \
			{catalog_path}""", step_number=3)

		print(f"Stats for {catalog_path}")
		run(f"python3 -m str_analysis.compute_catalog_stats --verbose {filtered_catalog_path}", step_number=3)

	# NOTE: we don't use the --merge-adjacent-loci-with-same-motif  option for str_analysis.merge_loci because
	# it's important to presenve locus definitions as they appear in the individual source catalogs. If loci
	# are merged together, this can create incompatibility with loci in source catalogs (ie. the illumina catalog) or
	# between future versions of the overall Simple Repeat Catalog since newly-added loci could merge with previously
	# added loci, causing loss of those loci from new vesions of the catalog. Variation cluster analysis is a better
	# way to merge adjacent loci where needed.
	catalog_paths = " ".join([
		f"{catalog_name}:{filtered_source_catalog_paths[catalog_name]}" for catalog_name, _ in source_catalogs_in_order
	])

	run(f"""python3 -u -m str_analysis.merge_loci --verbose \
		--add-found-in-fields \
		--output-format JSON \
		--discard-extra-fields-from-input-catalogs \
		--overlapping-loci-action keep-first \
		--write-merge-stats-tsv \
		--write-outer-join-table \
		--write-bed-files-with-unique-loci \
		--outer-join-overlap-table-min-sources 1 \
		--output-prefix {output_prefix}.merged \
		{catalog_paths}""", step_number=5)

	annotated_catalog_path = f"{output_prefix}.EH.with_annotations.json.gz"
	run(f"""python3 -u -m str_analysis.annotate_and_filter_str_catalog --verbose --show-progress-bar \
		--reference-fasta {args.hg38_reference_fasta} \
		--gene-models-source gencode \
		--gene-models-source refseq \
		--gene-models-source mane \
		--known-disease-associated-loci {primary_disease_associated_loci_path} \
		--min-motif-size {min_motif_size} \
		--max-motif-size {max_motif_size} \
		--min-interval-size-bp 1 \
		--discard-overlapping-intervals-with-similar-motifs \
		--output-path {annotated_catalog_path} \
		{output_prefix}.merged.json.gz""", step_number=6)

	# create a version of the ExpansionHunter catalog without extra annotations
	run(f"""python3 << EOF
import gzip, ijson, json

f = gzip.open("{annotated_catalog_path}", "rt")
out = gzip.open("{output_prefix}.EH.json.gz", "wt")
i = 0
out.write("[")
for record in ijson.items(f, "item", use_float=True):
	# skip chrM loci because ExpansionHunter prints an error like: 'Unable to extract chrM:-793-207 from hg38.fa'
	is_chrM = record["LocusId"].startswith("M-") or record["LocusId"].startswith("chrM-")
	if is_chrM: print(f"Skipping chrM locus: {{record['LocusId']}}")
	if is_chrM: continue
	if i > 0: out.write(", ")
	i += 1
	out.write(json.dumps({{ 
		k: v for k, v in record.items() if k in {{"LocusId", "ReferenceRegion", "VariantType", "LocusStructure"}} 
	}}, indent=4))
out.write("]")
EOF
""", step_number=7)

	run(f"python3 -m str_analysis.filter_out_loci_with_Ns_in_flanks "
		f"-R {args.hg38_reference_fasta} "
		f"-o {output_prefix}.EH.without_loci_with_Ns_in_flanks.json.gz "
		f"--output-list-of-filtered-loci {output_prefix}.loci_with_Ns_in_flanks.txt "
		f"{output_prefix}.EH.json.gz", step_number=8)
	run(f"mv {output_prefix}.EH.without_loci_with_Ns_in_flanks.json.gz {output_prefix}.EH.json.gz", step_number=8)

	release_files = [
		f"{output_prefix}.bed.gz",
		f"{output_prefix}.bed.gz.tbi",
		f"{annotated_catalog_path}",
		f"{output_prefix}.EH.json.gz",
		f"{output_prefix}.TRGT.bed",
		f"{output_prefix}.LongTR.bed",
		f"{output_prefix}.HipSTR.bed",
		f"{output_prefix}.GangSTR.bed",
	]


	# add variation cluster annotations to the catalog
	if args.variation_clusters_bed:
		variation_clusters_release_filename = f"{args.variation_clusters_output_prefix}.TRGT.bed.gz"
		variation_clusters_and_isolated_TRs_release_filename = args.variation_clusters_output_prefix.replace(
			"variation_clusters", "variation_clusters_and_isolated_TRs") + ".TRGT.bed.gz"

		assert variation_clusters_and_isolated_TRs_release_filename != variation_clusters_release_filename

		run(f"""python3 {base_dir}/scripts/add_variation_cluster_annotations_to_catalog.py \
			--verbose \
			--output-catalog-json-path {output_prefix}.EH.with_annotations.with_variation_clusters.json.gz \
			--known-pathogenic-loci-json-path {source_catalog_paths['KnownDiseaseAssociatedLoci']} \
			{args.variation_clusters_bed} \
			{annotated_catalog_path}""", step_number=9)

		run(f"mv {output_prefix}.EH.with_annotations.with_variation_clusters.json.gz {annotated_catalog_path}", step_number=9)
		run(f"cp {args.variation_clusters_bed} {variation_clusters_release_filename}", step_number=9)

		run(f"""python3 {base_dir}/scripts/add_isolated_loci_to_variation_cluster_catalog.py \
			--known-pathogenic-loci-json-path {source_catalog_paths['KnownDiseaseAssociatedLoci']} \
			-o {variation_clusters_and_isolated_TRs_release_filename} \
			{args.variation_clusters_bed} \
			{annotated_catalog_path}""", step_number=9)

		release_files.append(variation_clusters_release_filename)
		release_files.append(variation_clusters_and_isolated_TRs_release_filename)

		run(f"python3 {base_dir}/scripts/convert_trgt_catalog_to_longtr_format.py "
			f"--known-pathogenic-loci-json-path {source_catalog_paths['KnownDiseaseAssociatedLoci']} "
			f"{variation_clusters_release_filename}", step_number=10)
		run(f"python3 {base_dir}/scripts/convert_trgt_catalog_to_longtr_format.py "
			f"--known-pathogenic-loci-json-path {source_catalog_paths['KnownDiseaseAssociatedLoci']} "
			f"{variation_clusters_and_isolated_TRs_release_filename}", step_number=10)

		release_files.append(variation_clusters_release_filename.replace(".TRGT.bed.gz", ".LongTR.bed.gz"))
		release_files.append(variation_clusters_and_isolated_TRs_release_filename.replace(".TRGT.bed.gz", ".LongTR.bed.gz"))

	# add allele frequencies to the catalog
	run(f"""python3 -u {base_dir}/scripts/add_allele_frequency_annotations.py \
			--add-t2t-assembly-frequencies-to-overlapping-loci \
			-o {annotated_catalog_path}.with_allele_frequencies.json.gz  {annotated_catalog_path}""", step_number=11)

	run(f"mv {annotated_catalog_path}.with_allele_frequencies.json.gz {annotated_catalog_path}", step_number=11)

	# add LPS annotations
	if args.lps_annotations:
		run(f"""python3 {base_dir}/scripts/add_LPS_stdev_annotations_to_catalog.py \
			--output-catalog-json-path {output_prefix}.EH.with_annotations.with_LPS_annotations.json.gz \
			--known-pathogenic-loci-json-path {source_catalog_paths['KnownDiseaseAssociatedLoci']} \
			{args.lps_annotations} \
			{annotated_catalog_path}""", step_number=12)

		run(f"mv {output_prefix}.EH.with_annotations.with_LPS_annotations.json.gz {annotated_catalog_path}", step_number=12)

	# annotate with "TRsInRegion" based on adjacent loci
	if motif_size_label == "1_to_1000bp_motifs":
		adjacent_repeats_source_bed = f"{output_prefix}.bed.gz"

	# convert to BED
	run(f"python3 -m str_analysis.convert_expansion_hunter_catalog_to_bed --split-adjacent-repeats "
		f"{annotated_catalog_path}  --output-file {output_prefix}.bed.gz", step_number=13)

	run(f"python3 -m str_analysis.add_adjacent_loci_to_expansion_hunter_catalog "
		f"--ref-fasta {args.hg38_reference_fasta} "
		f"--source-of-adjacent-loci {adjacent_repeats_source_bed} "
		f"--add-extra-field TRsInRegion "
		f"--only-add-extra-fields "
		f"-o {annotated_catalog_path}.with_adjacent_loci_annotation.json.gz "
		f"{annotated_catalog_path}", step_number=14)

	run(f"mv {annotated_catalog_path}.with_adjacent_loci_annotation.json.gz {annotated_catalog_path}", step_number=14)


	# convert to TSV
	run(f"""python3 << EOF
import gzip
import json
import pandas as pd
from pprint import pformat
with gzip.open("{annotated_catalog_path}", "rt") as f:
	data = json.load(f)
print(f"Loaded {{len(data):,d}} records from {annotated_catalog_path}")
print("Writing it to TSV file...")
df = pd.DataFrame(data)
core_columns = [
	'LocusId', 'ReferenceRegion', 'LocusStructure', 'CanonicalMotif', 'TRsInRegion',
	'Source', 'GencodeGeneRegion', 'GencodeGeneId', 'GencodeGeneName', 'GencodeTranscriptId',
	'RefseqGeneRegion', 'RefseqGeneId', 'RefseqGeneName', 'RefseqTranscriptId', 
	'ManeGeneRegion', 'ManeGeneId', 'ManeGeneName', 'ManeTranscriptId',
	'KnownDiseaseAssociatedMotif',  'KnownDiseaseAssociatedLocus', 'NsInFlanks',
	'LeftFlankMappability', 'FlanksAndLocusMappability', 'RightFlankMappability', 
	'FoundInKnownDiseaseAssociatedLoci', 'FoundInIllumina174kPolymorphicTRs', 
	'FoundInPerfectRepeatsInReference', 'FoundInPolymorphicTRsInT2TAssemblies', 
	'NumRepeatsInReference', 'ReferenceRepeatPurity', 
	'AlleleFrequenciesFromIllumina174k', 'StdevFromIllumina174k',
	'AlleleFrequenciesFromT2TAssemblies', 'StdevFromT2TAssemblies',
	'VariationCluster', 'VariationClusterSizeDiff',
	'LPSLengthStdevFromHPRC100', 'LPSMotifFractionFromHPRC100',
]

drop_columns = ['VariantType', ]
for c in set(core_columns)  - set(df.columns): df[c] = None
df = df[core_columns + [c for c in df.columns if c not in (core_columns + drop_columns)]]

output_tsv_path = "{annotated_catalog_path.replace('.json.gz', '') + '.tsv.gz'}"
df.to_csv(output_tsv_path, sep="\\t", index=False)
print(f"Wrote {{len(df):,d}} rows to {{output_tsv_path}} with columns: {{pformat(list(df.columns))}}")
EOF
""", step_number=15)

	# convert the catalog from ExpansionHunter catalog format to TRGT, LongTR, HipSTR, and GangSTR formats
	run(f"python3 -m str_analysis.convert_expansion_hunter_catalog_to_trgt_catalog --split-adjacent-repeats {annotated_catalog_path}  --output-file {output_prefix}.TRGT.bed", step_number=16)
	run(f"python3 -m str_analysis.convert_expansion_hunter_catalog_to_longtr_format  {annotated_catalog_path}  --output-file {output_prefix}.LongTR.bed", step_number=17)
	run(f"python3 -m str_analysis.convert_expansion_hunter_catalog_to_hipstr_format  {annotated_catalog_path}  --output-file {output_prefix}.HipSTR.bed", step_number=18)
	run(f"python3 -m str_analysis.convert_expansion_hunter_catalog_to_gangstr_spec   {annotated_catalog_path}  --output-file {output_prefix}.GangSTR.bed", step_number=19)

	# Confirm that the TRGT catalog passes 'trgt validate'
	run(f"trgt validate --genome {args.hg38_reference_fasta}  --repeats {output_prefix}.TRGT.bed", step_number=20)

	# Perform basic internal consistency checks on the JSON catalog
	run(f"python3 {base_dir}/scripts/validate_catalog.py " +
		f"--known-pathogenic-loci-json-path {source_catalog_paths['KnownDiseaseAssociatedLoci']} " +
		("--check-for-presence-of-annotations --check-for-presence-of-all-known-loci " if motif_size_label == "1_to_1000bp_motifs" else "") +
		f"{annotated_catalog_path}", step_number=21)

	# copy files to the release_draft folder and compute catalog stats
	updated_release_files = []
	for path in release_files:
		if path.endswith(".bed"):
			if not os.path.isfile(f"{path}.gz"):
				if path.endswith(".TRGT.bed"):
					run(f"gzip -f {path}", step_number=22)  # TRGT v1.1.1 and lower only works with gzip, not bgzip
				else:
					run(f"bgzip -f {path}", step_number=22)
			updated_release_files.append(f"{path}.gz")
		else:
			if path.endswith(".json") or path.endswith(".json.gz") and ".EH." in path:
				run(f"python3 {base_dir}/scripts/validate_json.py -k LocusId -k LocusStructure -k ReferenceRegion -k VariantType {path}", step_number=22)
			updated_release_files.append(path)

	if release_tar_gz_path is None:
		for path in updated_release_files:
			run(f"cp {path} {release_draft_folder}", step_number=22)
	else:
		run(f"tar czf {release_tar_gz_path} -C {os.path.dirname(output_prefix)} " + " ".join([os.path.basename(p) for p in updated_release_files]), step_number=22)
		run(f"cp {release_tar_gz_path} {release_draft_folder}", step_number=22)

	run(f"python3 -m str_analysis.compute_catalog_stats --verbose {annotated_catalog_path}", step_number=23)

	# report hours, minutes, seconds relative to start_time
	diff = time.time() - start_time
	print(f"Done generating {output_prefix} catalog. Took {diff//3600:.0f}h, {(diff%3600)//60:.0f}m, {diff%60:.0f}s")

	if motif_size_label != "1_to_1000bp_motifs":
		continue

	# compare to the GangSTR_v17 catalog to make sure it's included
	comparison_catalogs_in_order = [
		("GangSTR_v17", "https://s3.amazonaws.com/gangstr/hg38/genomewide/hg38_ver17.bed.gz"),
		#("vamos_catalog_v2.1", "https://storage.googleapis.com/str-truth-set/hg38/ref/other/vamos_catalog.v2.1.bed.gz"),
	]

	comparison_catalog_paths = {}
	for catalog_name, url in comparison_catalogs_in_order:
		if not os.path.isfile(os.path.basename(url)):
			run(f"wget -O {os.path.basename(url)}.tmp -qnc {url} && mv {os.path.basename(url)}.tmp {os.path.basename(url)}")
		comparison_catalog_paths[catalog_name] = os.path.abspath(os.path.basename(url))

	path_after_conversion = comparison_catalog_paths["GangSTR_v17"].replace(".bed.gz", ".json.gz")
	run(f"python3 -u -m str_analysis.convert_gangstr_spec_to_expansion_hunter_catalog --verbose {comparison_catalog_paths['GangSTR_v17']} -o {path_after_conversion}", step_number=30)
	comparison_catalog_paths["GangSTR_v17"] = path_after_conversion

	# compare catalog to other catalogs
	start_time = time.time()

	for catalog_name, path in comparison_catalog_paths.items():
		filtered_comparison_catalog_path = re.sub("(.json|.bed)(.gz)?$", "", path) + ".filtered.json.gz"

		run(f"""python3 -m str_analysis.annotate_and_filter_str_catalog \
			--reference-fasta {args.hg38_reference_fasta} \
			--skip-gene-annotations \
			--skip-mappability-annotations \
			--skip-disease-loci-annotations \
			--min-motif-size {min_motif_size} \
			--max-motif-size {max_motif_size} \
			--output-path {filtered_comparison_catalog_path} \
			--verbose \
			{path}""", step_number=31)

		run(f"python3 -m str_analysis.compute_catalog_stats --verbose {filtered_comparison_catalog_path}", step_number=32)

		run(f"""python3 -u -m str_analysis.merge_loci \
			--output-prefix {catalog_name} \
			--output-format JSON \
			--overlapping-loci-action keep-first \
			--verbose \
			--write-merge-stats-tsv \
			{annotated_catalog_path} \
			{filtered_comparison_catalog_path}""", step_number=33)

	diff = time.time() - start_time
	print(f"Done with comparisons. Took {diff//3600:.0f}h, {(diff%3600)//60:.0f}m, {diff%60:.0f}s")

