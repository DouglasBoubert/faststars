#! /bin/bash
set -e

PATH=/opt/local/bin:/usr/local/bin:$PATH ; export PATH
LD_LIBRARY_PATH=/usr/local/lib:/opt/local/lib ; export LD_LIBRARY_PATH

cd /var/www/html/hvs/astrocats
python -m astrocats.scripts.webcat -c hvs &
pids[0]=$!
python -m astrocats.scripts.webcat -c hvs -by &
pids[1]=$!
python -m astrocats.faststars.scripts.dupecat &
pids[2]=$!
python -m astrocats.faststars.scripts.conflictcat &
pids[3]=$!
python -m astrocats.scripts.bibliocat -c hvs &
pids[4]=$!
python -m astrocats.faststars.scripts.erratacat &
pids[5]=$!
#python -m astrocats.scripts.hostcat -c hvs &
#pids[6]=$!
python -m astrocats.scripts.hammertime -c hvs &
pids[7]=$!
#python -m astrocats.faststars.scripts.histograms &
#pids[8]=$!
#python -m astrocats.scripts.atelscbetsiaucs -c hvs &
#pids[9]=$!
#python -m astrocats.faststars.scripts.frbcat &
#pids[10]=$!
for pid in ${pids[*]}; do
	wait $pid
done
cd /var/www/html/hvs/astrocats/astrocats/faststars/output/html
bash thumbs.sh
cd /var/www/html/hvs/astrocats
