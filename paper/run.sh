set -ex

nohup python3 -u compare_catalogs.py --hg38 ~/hg38.fa --force
cat nohup.out | grep -v collaps > compare_catalogs.log

#nohup python3 -u compare_catalogs.py --hg38 ~/hg38.fa
#cat nohup.out | grep -v collaps > compare_catalogs.stats_only.log
