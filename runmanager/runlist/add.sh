#!/bin/sh

for i in $(seq 6000 8000); do
#for i in $(seq 7164 8000); do
#for i in $(seq 6200 7000); do
#for i in $(seq 6106 7000); do
# for i in $(seq 5284 5900); do
    r="/group/had/sks/E40/JPARC2019Feb/e40_2019feb/run0${i}.dat.gz"
    #r=`find /group/had/sks/E40/JPARC2018Jun/e40_2018jun/ -name "*${i}*.dat.gz"`
    if [ -f $r ]; then
	#echo "  $i:" >> mstmonitor_all.yml
	echo "  $i:" >> scaler_2019feb.yml
	#echo "  $i:" >> scaler_2018jun.yml
    fi
done
