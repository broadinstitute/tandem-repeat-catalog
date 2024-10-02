rm nohup.out

set -ex

nohup python3 -u compare_catalogs.py --hg38 ~/hg38.fa 
cp nohup.out compare_catalogs.log

