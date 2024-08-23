set -ex

#nohup python3 -u compare_catalogs.py --hg38 ~/hg38.fa --force
#cat nohup.out  > compare_catalogs.log

nohup python3 -u compare_catalogs.py --outer-join --hg38 ~/hg38.fa
cat nohup.out  > compare_catalogs.stats_only.log
