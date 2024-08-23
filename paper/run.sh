rm nohup.out

set -ex

nohup python3 -u compare_catalogs.py --hg38 ~/hg38.fa  --outer-join --force-annotate  
cat nohup.out  > compare_catalogs.log

#nohup python3 -u compare_catalogs.py --hg38 ~/hg38.fa  --outer-join --force-stats
#cat nohup.out  > compare_catalogs.stats_only.log

#nohup python3 -u compare_catalogs.py --hg38 ~/hg38.fa  --outer-join
#cat nohup.out  > compare_catalogs.stats_only.log
