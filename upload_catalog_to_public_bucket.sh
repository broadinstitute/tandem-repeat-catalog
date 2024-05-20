set -ex

gsutil -m cp README.md merged_and_annotated_catalog.hg38.* gs://str-truth-set/hg38/ref/other/colab-repeat-finder/hg38_repeats.motifs_1_to_50bp.repeats_3x_and_spans_9bp/
