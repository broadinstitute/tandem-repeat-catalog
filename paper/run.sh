set -ex

nohup python3 -u compare_catalogs.py --hg38 ~/hg38.fa
#cat nohup.out | grep -v collaps > all_steps.merge_and_annotate_loci.log
