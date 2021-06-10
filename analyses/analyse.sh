#!/bin/bash
cd $1
mkdir -p output/img/ForagingModel
mkdir -p output/mp4
echo "running node"
time node ../../src/runForagingModel-$1chemotaxis.js ../../results/results-$1-$2/evolutionary-potts/params/gen30/$3.json true > $2_$3.csv
echo "running done"
head -n -1 $2_$3.csv > $2_$3_temp.csv; mv $2_$3_temp.csv $2_$3.csv
mv output $2_$3_o
../../makemovie.sh $2_$3_o