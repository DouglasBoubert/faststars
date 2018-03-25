#! /bin/bash
set -e

PATH=/root/miniconda3/bin:/opt/local/bin:/usr/local/bin:$PATH ; export PATH
LD_LIBRARY_PATH=/usr/local/lib:/opt/local/lib ; export LD_LIBRARY_PATH

cd /var/www/html/hvs/astrocats
python -m astrocats faststars git-pull
python -m astrocats faststars import
SNEUPDATE=$?
echo $SNEUPDATE
if [[ $SNEUPDATE == 0 ]]; then
	astrocats/faststars/scripts/generate-web.sh
	python -m astrocats faststars git-pull
	python -m astrocats faststars git-push
	#stamp=$(date +"%Y-%m-%d %k:%M")
	#./commit-and-push-repos.sh "Auto update: $stamp"
fi
