#!/bin/bash
data_type=$1
MAXEXPNUM=$2
K=$3
python Drawing/find_best_K.py ${data_type} ${MAXEXPNUM} ${K}
FILENAME=buf/${data_type}.txt
cnt=0
array=()
while read line;
do
    cnt=$(expr $cnt + 1)
    if test $cnt -eq 1; then
	num_data=$line
    fi
    if test $cnt -eq 2; then
	number_of_topic=$line
    fi
    if test $cnt -ge 3; then
	array+=($line)
    fi
done<$FILENAME
for i in $(seq 2 ${K}); do
    d_num=${array[i-2]}
    str1=result/data${data_type}_${d_num}
    str2=../data${data_type}_${num_data}
    cd ${str1}
    if test ${i} -le 9; then
	str3=result_k0${i}.txt
    else
	str3=result_k${i}.txt
    fi
    cp ${str3} ${str2}
    cd ../..
done
for i in $(seq 1 ${MAXEXPNUM}); do
    str1=result/data${data_type}_${i}
    if test ${i} -eq ${num_data}; then
	cd ${str1}
	mkdir figure
	cd ../..
	str2=result/data${data_type}
	mv ${str1} ${str2}
    else
	rm -r ${str1}
    fi
done
python Drawing/comparison_K.py ${data_type} ${K}
str2=result/data${data_type}
cd ${str2}
cd figure
mkdir ${number_of_topic}_signature
mkdir ${number_of_topic}_arrangement
cd ../../..
bash scripts/draw_signature.sh ${data_type} ${number_of_topic}
