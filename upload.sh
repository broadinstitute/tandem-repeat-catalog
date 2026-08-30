set -ex

gcloud storage cp results__2026-08-30/release_draft_2026-08-30/*.*   gs://tandem-repeat-catalog/v2.1/

# The extended locus boundaries source catalog is generated separately by
# extended_loci/generate_extended_loci_catalog.sh and is the default input for the pipeline's
# --extended-loci-catalog argument, so it is published alongside the release it belongs to.
gcloud storage cp extended_loci/TRExplorer.extended_loci_v2.1.hg38.bed.gz*   gs://tandem-repeat-catalog/v2.1/
