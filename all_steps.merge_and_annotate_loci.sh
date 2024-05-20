set -ex

REFERENCE_FASTA_PATH=${1:-~/hg38.fa}

# check if REFERENCE_FASTA_PATH exists
if [ ! -f ${REFERENCE_FASTA_PATH} ]; then
  echo "ERROR: GRCh38 reference fasta file not found @ ${REFERENCE_FASTA_PATH}. Please specify the correct path as a command-line arg."
  exit 1
fi

# STEP #1  (already done for hg38, and results are publicly available, so just download them in step #2)

echo Skipping step1. Will just download the results from the public URL.
#git clone git@github.com:broadinstitute/colab-repeat-finder.git
#cd colab-repeat-finder/python
#python3 perfect_repeat_finder.py  --min-repeats 3  --min-span 9  --min-motif-size 1  --max-motif-size 50  --output-prefix perfect_repeats.hg38  --show-progress-bar  ${REFERENCE_FASTA_PATH}


# STEP #2

KNOWN_DISEASE_ASSOCIATED_LOCI=https://raw.githubusercontent.com/broadinstitute/str-analysis/main/str_analysis/variant_catalogs/variant_catalog_without_offtargets.GRCh38.json
ALL_PERFECT_REPEATS_CATALOG_URL=https://storage.googleapis.com/str-truth-set/hg38/ref/other/colab-repeat-finder/hg38_repeats.motifs_1_to_50bp.repeats_3x_and_spans_9bp/hg38_repeats.motifs_1_to_50bp.repeats_3x_and_spans_9bp.bed.gz
ILLUMINA_CATALOG_URL=https://storage.googleapis.com/str-truth-set/hg38/ref/other/illumina_variant_catalog.sorted.bed.gz
TRUTH_SET_CATALOG_URL=https://storage.googleapis.com/str-truth-set-v2/filter_vcf/all_repeats_including_homopolymers_keeping_loci_that_have_overlapping_variants/combined/combined.51_samples.positive_loci.json

wget -nc $KNOWN_DISEASE_ASSOCIATED_LOCI
wget -nc $ALL_PERFECT_REPEATS_CATALOG_URL
wget -nc $ILLUMINA_CATALOG_URL
wget -nc $TRUTH_SET_CATALOG_URL

time python3 -u -m str_analysis.merge_loci --verbose \
  --merge-adjacent-loci-with-same-motif \
  --add-extra-fields-from-input-catalogs \
  --output-format JSON \
  --output-prefix merged_catalog \
  $(basename ${KNOWN_DISEASE_ASSOCIATED_LOCI}) \
  $(basename ${ALL_PERFECT_REPEATS_CATALOG_URL}) \
  $(basename ${ILLUMINA_CATALOG_URL}) \
  $(basename ${TRUTH_SET_CATALOG_URL})

time python3 -u -m str_analysis.annotate_and_filter_str_catalog --verbose \
   --skip-gene-annotations \
   --skip-disease-loci-annotations \
   --skip-mappability-annotations \
   --discard-loci-with-non-acgt-bases \
   --reference ${REFERENCE_FASTA_PATH} \
   --output-path merged_and_annotated_catalog.hg38.json.gz \
   merged_catalog.json.gz


# STEP #3

python3 -m str_analysis.convert_expansion_hunter_variant_catalog_to_trgt_catalog   merged_and_annotated_catalog.json.gz  --output-file merged_and_annotated_catalog.hg38.TRGT.bed
python3 -m str_analysis.convert_expansion_hunter_variant_catalog_to_longtr_format  merged_and_annotated_catalog.json.gz  --output-file merged_and_annotated_catalog.hg38.LongTR.bed
python3 -m str_analysis.convert_expansion_hunter_variant_catalog_to_hipstr_format  merged_and_annotated_catalog.json.gz  --output-file merged_and_annotated_catalog.hg38.HipSTR.bed
python3 -m str_analysis.convert_expansion_hunter_variant_catalog_to_gangstr_spec   merged_and_annotated_catalog.json.gz  --output-file merged_and_annotated_catalog.hg38.GangSTR.bed
