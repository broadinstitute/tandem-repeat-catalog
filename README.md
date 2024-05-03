### Authors

* Ben Weisburd, Broad Institute
* Egor Dolzhenko, PacBio

### Introduction

Tandem repeats are regions of the genome that consist of consecutive copies of some motif sequence. For example, CAGCAGCAG is a tandem repeat with motif CAG. Many types of genomic studies require annotations of tandem repeats in the reference genome, called repeat catalogs. Repeat catalogs typically consist of the repeat reference coordinates and one or multiple motif sequences that the repeat is composed of.

The purpose of this repo is to provide methods and best practices for defining high-quality repeat catalogs. Although our initial focus is on the human genome, this tutorial should be also applicable to other closely-related genomes. We'd also love to extend this work to plants and other species. Please consider creating a GitHub issue or reaching out by email if you are interested in this.

### Defining a genome-wide tandem repeat catalog

The steps below start with a reference genome fasta file and proceed to generate a comprehensive genome-wide TR catalog that can be used with the following TR genotyping tools:
* [ExpansionHunter](https://github.com/Illumina/ExpansionHunter)
* [TRGT](https://github.com/PacificBiosciences/trgt)
* [LongTR](https://github.com/gymrek-lab/LongTR)
* [HipSTR](https://github.com/HipSTR-Tool/HipSTR)
* [GangSTR](https://github.com/gymreklab/GangSTR)


### Step 1: Detect all perfect (ie. non-interrupted) tandem repeats in the reference genome

Run [colab-repeat-finder](https://github.com/broadinstitute/colab-repeat-finder) to detect perfect (ie. non-interrupted) tandem repeats in the reference genome that satisfy the following two criteria:
* Have a repeat motif that's between 1bp and 50bp, including all STRs (1-6bp) and many VNTR sizes (7-50bp)
* Have 3 or more consecutive repeats in the reference genome that span at least 9bp 

The output for hg38 is already available [here](https://storage.googleapis.com/str-truth-set/hg38/ref/other/colab-repeat-finder/hg38_repeats.motifs_1_to_50bp.repeats_3x_and_spans_9bp/hg38_repeats.motifs_1_to_50bp.repeats_3x_and_spans_9bp.bed.gz).

**Commands:**
```
REFERENCE_FASTA_PATH=/path/to/hg38.fa

git clone git@github.com:broadinstitute/colab-repeat-finder.git
cd colab-repeat-finder/python
python3 perfect_repeat_finder.py  --min-repeats 3  --min-span 9  --min-motif-size 1  --max-motif-size 50  --output-prefix perfect_repeats.hg38  --show-progress-bar  ${REFERENCE_FASTA_PATH}
```

*NOTE:* colab-repeat-finder is a new, simple tool that finds all perfect repeats in a given input sequence that satify user parameters. If you are interested in profiling imperfect repeats, you can generate the repeat catalog using TandemRepeatFinder (TRF).


### Step 2: Merge and annotate loci from step 1 (and optionally add other repeat catalogs)

Here, we fix some representation issues in the catalog from step 1:
- collapse adjacent loci that have the same motif but were reported as two separate loci by colab-repeat-finder due to a single base pair interruption between them.
- filter out repeats that include non-A,C,G,T characters, including Ns

Optionally, we can augment our catalog with loci from other sources. Specifically, while the catalog from step 1 includes all perfect tandem repeats loci that have at least 3 repeats of some motif and span at least 9bp in the reference genome, it misses polymorphic loci that coincidentally have 2 or fewer repeats in the reference (with larger TR alleles segregating in the population). To capture these loci, we can merge the catalog from step 1 with any available catalogs of polymorphic tandem repeat loci that were generated via orthogonal methods. For the human genome, these include:

* [Illumina catalog](https://github.com/Illumina/RepeatCatalogs) of 174k TR loci that are polymporphic in 2.4k diverse population samples from 1kGP.  
* Truth set from [[weisburd 2023](https://www.biorxiv.org/content/10.1101/2023.05.05.539588v1)] that has now been updated to include polymorphic TR loci in 51 samples from the HPRC. 


The following commands should be run even if you only have the one catalog from step1 and did not add any other catalogs. To do this, simply remove the additional catalogs from the commmands below.

**Commands:**
```
ALL_PERFECT_REPEATS_CATALOG_URL=https://storage.googleapis.com/str-truth-set/hg38/ref/other/colab-repeat-finder/hg38_repeats.motifs_1_to_50bp.repeats_3x_and_spans_9bp/hg38_repeats.motifs_1_to_50bp.repeats_3x_and_spans_9bp.bed.gz
ILLUMINA_CATALOG_URL=https://storage.googleapis.com/str-truth-set/hg38/ref/other/illumina_variant_catalog.sorted.bed.gz
TRUTH_SET_CATALOG_URL=https://storage.googleapis.com/str-truth-set-v2/filter_vcf/all_repeats_including_homopolymers_keeping_loci_that_have_overlapping_variants/combined/combined.51_samples.positive_loci.json


wget $ALL_PERFECT_REPEATS_CATALOG_URL
wget $ILLUMINA_CATALOG_URL
wget $TRUTH_SET_CATALOG_URL

python3 -u -m str_analysis.merge_loci --verbose \
  --merge-adjacent-loci-with-same-motif \
  --output-format JSON \
  --output-prefix merged_catalog \
  $(basename ${ALL_PERFECT_REPEATS_CATALOG_URL}) \
  $(basename ${ILLUMINA_CATALOG_URL}) \
  $(basename ${TRUTH_SET_CATALOG_URL})

python3 -u -m str_analysis.annotate_and_filter_str_catalog --verbose \
   --skip-gene-annotations \
   --skip-disease-loci-annotations \
   --skip-mappability-annotations \
   --discard-loci-with-non-acgt-bases \
   --reference ${REFERENCE_FASTA_PATH} \
   --output-path merged_and_annotated_catalog.json.gz \
   merged_catalog.json.gz
```

### Step 3: Convert the combined catalog into tool-specific catalog formats

The `merged_catalog.json.gz` file from step2 is already in a format that can be gunzipped and passed as input to [ExpansionHunter](https://github.com/Illumina/ExpansionHunter) (or to [this fork of ExpansionHunter](https://github.com/bw2/ExpansionHunter) which can read the gzipped JSON file directly). The commands below convert `merged_catalog.json.gz` to formats which can be passed as input to:
* [TRGT](https://github.com/PacificBiosciences/trgt)
* [LongTR](https://github.com/gymrek-lab/LongTR)
* [HipSTR](https://github.com/HipSTR-Tool/HipSTR)
* [GangSTR](https://github.com/gymreklab/GangSTR)

**Commands:**

```
python3 -m str_analysis.convert_expansion_hunter_variant_catalog_to_trgt_catalog   --output-file merged_and_annotated_catalog.TRGT.bed     merged_and_annotated_catalog.json.gz  
python3 -m str_analysis.convert_expansion_hunter_variant_catalog_to_longtr_format  --output-file merged_and_annotated_catalog.LongTR.bed   merged_and_annotated_catalog.json.gz  
python3 -m str_analysis.convert_expansion_hunter_variant_catalog_to_hipstr_format  --output-file merged_and_annotated_catalog.HipSTR.bed   merged_and_annotated_catalog.json.gz  
python3 -m str_analysis.convert_expansion_hunter_variant_catalog_to_gangstr_spec   --output-file merged_and_annotated_catalog.GangSTR.bed  merged_and_annotated_catalog.json.gz  
```

