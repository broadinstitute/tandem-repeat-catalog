set -ex

# install the latest version of the str-analysis package
#python3 -m pip install --upgrade git+https://github.com/broadinstitute/str-analysis.git

REFERENCE_FASTA_PATH=~/hg38.fa
GENCODE_GTF_PATH=gencode.v46.basic.annotation.gtf.gz
# STEP #1  (already done for hg38, and results are publicly available. Will just download them in step #2)

#git clone git@github.com:broadinstitute/colab-repeat-finder.git
#cd colab-repeat-finder/python
#python3 perfect_repeat_finder.py  --min-repeats 3  --min-span 9  --min-motif-size 1  --max-motif-size 50  --output-prefix perfect_repeats.hg38  --show-progress-bar  ${REFERENCE_FASTA_PATH}


# STEP #2:  generate catalog
KNOWN_DISEASE_ASSOCIATED_LOCI_URL=https://raw.githubusercontent.com/broadinstitute/str-analysis/main/str_analysis/variant_catalogs/variant_catalog_without_offtargets.GRCh38.json
#ALL_PERFECT_REPEATS_CATALOG_URL=https://storage.googleapis.com/str-truth-set/hg38/ref/other/colab-repeat-finder/hg38_repeats.motifs_1_to_50bp.repeats_3x_and_spans_9bp/hg38_repeats.motifs_1_to_50bp.repeats_3x_and_spans_9bp.bed.gz
ALL_PERFECT_REPEATS_CATALOG_URL=https://storage.googleapis.com/str-truth-set/hg38/ref/other/colab-repeat-finder/hg38_repeats.motifs_1_to_1000bp.repeats_3x_and_spans_9bp/hg38_repeats.motifs_1_to_1000bp.repeats_3x_and_spans_9bp.bed.gz
ILLUMINA_CATALOG_URL=https://storage.googleapis.com/str-truth-set/hg38/ref/other/illumina_variant_catalog.sorted.bed.gz
TRUTH_SET_CATALOG_URL=https://storage.googleapis.com/str-truth-set-v2/filter_vcf/all_repeats_including_homopolymers_keeping_loci_that_have_overlapping_variants/combined/combined.51_samples.positive_loci.json

wget -nc $KNOWN_DISEASE_ASSOCIATED_LOCI_URL
wget -nc $ALL_PERFECT_REPEATS_CATALOG_URL
wget -nc $ILLUMINA_CATALOG_URL
wget -nc $TRUTH_SET_CATALOG_URL

KNOWN_DISEASE_ASSOCIATED_LOCI_PATH=$(basename ${KNOWN_DISEASE_ASSOCIATED_LOCI_URL})
ALL_PERFECT_REPEATS_CATALOG_PATH=$(basename ${ALL_PERFECT_REPEATS_CATALOG_URL})
ILLUMINA_CATALOG_PATH=$(basename ${ILLUMINA_CATALOG_URL})
TRUTH_SET_CATALOG_PATH=$(basename ${TRUTH_SET_CATALOG_URL})

python3 -m str_analysis.split_adjacent_loci_in_expansion_hunter_catalog ${KNOWN_DISEASE_ASSOCIATED_LOCI_PATH}

KNOWN_DISEASE_ASSOCIATED_LOCI_PATH=$(echo $KNOWN_DISEASE_ASSOCIATED_LOCI_PATH | sed 's/.json/.split.json/')

# change motif definition for the RFC1 locus from AARRG => AAAAG since the catalog doesn't currently support IUPAC codes
sed -i 's/AARRG/AAAAG/g' ${KNOWN_DISEASE_ASSOCIATED_LOCI_PATH}

set +x
for catalog_path in \
  ${KNOWN_DISEASE_ASSOCIATED_LOCI_PATH} \
  ${ALL_PERFECT_REPEATS_CATALOG_PATH} \
  ${ILLUMINA_CATALOG_PATH} \
  ${TRUTH_SET_CATALOG_PATH}
do

    filtered_catalog_path=$(echo $catalog_path | sed 's/.json$//' | sed 's/.bed.gz$//').filtered.json
    set -x
    python3 -u -m str_analysis.annotate_and_filter_str_catalog --verbose \
      --reference-fasta ${REFERENCE_FASTA_PATH} \
      --min-motif-size 2 \
      --max-motif-size 1000 \
      --min-interval-size-bp 1 \
      --skip-gene-annotations \
      --skip-mappability-annotations \
      --skip-disease-loci-annotations \
      --discard-loci-with-non-ACGT-bases-in-reference \
      --discard-loci-with-non-ACGTN-bases-in-motif \
      --output-path ${filtered_catalog_path} \
      ${catalog_path}

    set +x
    echo "Stats for ${catalog_path}"
    set -x
    python3 -m str_analysis.compute_catalog_stats --verbose ${filtered_catalog_path}
    set +x
done
set -x

KNOWN_DISEASE_ASSOCIATED_LOCI_FILTERED_PATH=$(echo $KNOWN_DISEASE_ASSOCIATED_LOCI_PATH | sed 's/.json/.filtered.json/')
ALL_PERFECT_REPEATS_CATALOG_FILTERED_PATH=$(echo $ALL_PERFECT_REPEATS_CATALOG_PATH | sed 's/.bed.gz/.filtered.json/')
ILLUMINA_CATALOG_FILTERED_PATH=$(echo $ILLUMINA_CATALOG_PATH | sed 's/.bed.gz/.filtered.json/')
TRUTH_SET_CATALOG_FILTERED_PATH=$(echo $TRUTH_SET_CATALOG_PATH | sed 's/.json/.filtered.json/')

python3 -u -m str_analysis.merge_loci --verbose \
	--merge-adjacent-loci-with-same-motif \
	--add-source-field \
	--output-format JSON \
	--write-bed-files-with-new-loci \
	--output-prefix repeat_catalog.3x_and_9bp.hg38.merged \
	${KNOWN_DISEASE_ASSOCIATED_LOCI_FILTERED_PATH} \
	${ALL_PERFECT_REPEATS_CATALOG_FILTERED_PATH} \
	${ILLUMINA_CATALOG_FILTERED_PATH} \
	${TRUTH_SET_CATALOG_FILTERED_PATH}

# debug merge
#for catalog_path in \
#  ${KNOWN_DISEASE_ASSOCIATED_LOCI_FILTERED_PATH} \
#  ${ALL_PERFECT_REPEATS_CATALOG_FILTERED_PATH} \
#  ${ILLUMINA_CATALOG_FILTERED_PATH} \
#  ${TRUTH_SET_CATALOG_FILTERED_PATH}
#do
#  echo =================================================
#  python3 -u -m str_analysis.merge_loci --verbose \
#    --write-bed-files-with-new-loci \
#    repeat_catalog.3x_and_9bp.hg38.merged.json.gz \
#    ${catalog_path}
#done


python3 -u -m str_analysis.annotate_and_filter_str_catalog --verbose \
	--reference-fasta ${REFERENCE_FASTA_PATH} \
	--genes-gtf ${GENCODE_GTF_PATH} \
	--gene-models-source gencode \
	--min-motif-size 2 \
	--max-motif-size 1000 \
	--min-interval-size-bp 1 \
	--discard-overlapping-intervals-with-similar-motifs \
	--output-path repeat_catalog.3x_and_9bp.hg38.merged_and_annotated.json.gz \
	repeat_catalog.3x_and_9bp.hg38.merged.json.gz


# STEP #3: convert the catalog from ExpansionHunter catalog format to TRGT, LongTR, HipSTR, and GangSTR formats
python3 -m str_analysis.convert_expansion_hunter_variant_catalog_to_trgt_catalog --split-adjacent-repeats repeat_catalog.3x_and_9bp.hg38.merged_and_annotated.json.gz  --output-file repeat_catalog.3x_and_9bp.hg38.TRGT.bed &
python3 -m str_analysis.convert_expansion_hunter_variant_catalog_to_longtr_format  repeat_catalog.3x_and_9bp.hg38.merged_and_annotated.json.gz  --output-file repeat_catalog.3x_and_9bp.hg38.LongTR.bed &
python3 -m str_analysis.convert_expansion_hunter_variant_catalog_to_hipstr_format  repeat_catalog.3x_and_9bp.hg38.merged_and_annotated.json.gz  --output-file repeat_catalog.3x_and_9bp.hg38.HipSTR.bed &
python3 -m str_analysis.convert_expansion_hunter_variant_catalog_to_gangstr_spec   repeat_catalog.3x_and_9bp.hg38.merged_and_annotated.json.gz  --output-file repeat_catalog.3x_and_9bp.hg38.GangSTR.bed &
wait




# STEP #4: create release_draft folder and compute catalog stats
timestamp=$(date '+%Y-%m-%d')
release_draft_folder=release_draft_${timestamp}
mkdir -p ${release_draft_folder}
for path in repeat_catalog.3x_and_9bp.hg38.merged_and_annotated.bed.gz \
	repeat_catalog.3x_and_9bp.hg38.merged_and_annotated.bed.gz.tbi \
	repeat_catalog.3x_and_9bp.hg38.merged_and_annotated.json.gz \
	repeat_catalog.3x_and_9bp.hg38.TRGT.bed.gz \
	repeat_catalog.3x_and_9bp.hg38.LongTR.bed \
	repeat_catalog.3x_and_9bp.hg38.HipSTR.bed \
	repeat_catalog.3x_and_9bp.hg38.GangSTR.bed;
do
    cp ${path} ${release_draft_folder}/
done

python3 -m str_analysis.compute_catalog_stats --verbose repeat_catalog.3x_and_9bp.hg38.merged_and_annotated.json.gz


# STEP #5: compare it to other catalogs
wget -nc https://s3.amazonaws.com/gangstr/hg38/genomewide/hg38_ver17.bed.gz
wget -nc https://storage.googleapis.com/str-truth-set/hg38/ref/other/adotto_tr_catalog_v1.2.bed.gz
wget -nc https://storage.googleapis.com/str-truth-set/hg38/ref/other/mukamel_VNTR_catalog.bed.gz
wget -nc https://storage.googleapis.com/str-truth-set/hg38/ref/other/vamos_catalog.v2.1.bed.gz

python3 -u -m str_analysis.convert_gangstr_spec_to_expansion_hunter_variant_catalog --verbose hg38_ver17.bed.gz

for catalog_filename in mukamel_VNTR_catalog.bed.gz  vamos_catalog.v2.1.bed.gz hg38_ver17.variant_catalog.json    adotto_tr_catalog_v1.2.bed.gz
do
    echo Processing ${catalog_filename}
    output_prefix=$(echo ${catalog_filename} | sed 's/.bed.gz$//' | sed 's/.variant_catalog.json$//')

    python3 -m str_analysis.annotate_and_filter_str_catalog \
	    -r ${REFERENCE_FASTA_PATH} \
	    --skip-gene-annotations \
	    --skip-mappability-annotations \
	    --skip-disease-loci-annotations \
	    --min-motif-size 2 \
	    --max-motif-size 1000 \
	    --output-path ${output_prefix}.annotated_and_filtered.json \
	    --verbose \
	    ${catalog_filename}

    python3 -m str_analysis.compute_catalog_stats --verbose ${output_prefix}.annotated_and_filtered.json

    python3 -u -m str_analysis.merge_loci --verbose  \
	    --output-prefix ${output_prefix}  \
	    --write-bed-files-with-new-loci \
	    repeat_catalog.3x_and_9bp.hg38.merged_and_annotated.json.gz \
	    ${output_prefix}.annotated_and_filtered.json
done
