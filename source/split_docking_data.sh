#!/bin/bash

root_dir=$1
splits_dir=$2
n_splits=$3

prot_count=0
declare qts=()
declare boxes=()
declare configs=()
for file in $(ls "$root_dir"/**/*.pdbqt); do
	qts+=($file)
	boxes+=(${file%.pdbqt}.box.pdb)
	configs+=(${file%.pdbqt}.box.txt)
	((prot_count++))
done

echo $prot_count
echo ${qts[0]}
echo ${boxes[0]}
echo ${configs[0]}

split_size=$(($prot_count / $n_splits))
current_split_size=0
split_idx=0
mkdir "$splits_dir"

out_dir="$splits_dir"/splits_"$split_idx"
mkdir "$out_dir"


for i in $(seq 0 $(($prot_count - 1))); do
	echo $i
    
	cp "${qts[i]}" "$out_dir"
	cp "${boxes[i]}" "$out_dir"
	cp "${configs[i]}" "$out_dir"
	((current_split_size++))
	if [ "$current_split_size" -eq "$split_size" ]; then
		((split_idx++))
		current_split_size=0
		out_dir="$splits_dir"/splits_"$split_idx"
		mkdir "$out_dir"
	fi
done
