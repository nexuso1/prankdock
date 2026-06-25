#!/bin/bash

root_dir=$1
splits_dir=$2
n_splits=$3

count=0
declare qts=()
declare boxes=()
declare configs=()
for file in $(ls "$root_dir"/**/*.pdbqt); do
        qts+=($file)
        boxes+=(${file%.pdbqt}.box.pdb)
        configs+=(${file%.pdbqt}.box.txt)
	count=$((count+1))
done

echo $count
echo ${qts[0]}
echo ${boxes[0]}
echo ${configs[0]}

size=$((count / $n_splits))

mkdir "$splits_dir"
for i in $(seq 0 $(($n_splits - 1))); do
	echo $i
    out_dir=./"$splits_dir"/splits_"$i"
	mkdir "$out_dir"
	for j in $(seq 0 $size); do
		index=$((i * size + j))
		cp "${qts[index]}" "$out_dir"
		cp "${boxes[index]}" "$out_dir"
		cp "${configs[index]}" "$out_dir"
	done
done

for i in $(seq $(($count - $count % $size)) $(($count - 1))); do
	echo $i
	s_idx=$((n_splits - 1))
	cp "${qts[i]}" "./splits_$s_idx/"
	cp "${boxes[i]}" "./splits_$s_idx/"
	cp "${configs[i]}" "./splits_$s_idx/"
done
