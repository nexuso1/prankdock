#!/bin/bash
count=0
declare pdbs=()
for file in $(ls ../data/pdbs/*.pdb); do
        pdbs+=($file)
	count=$((count+1))
done

echo $count
echo ${pdbs[0]}

n_splits=$1 
size=$((count / $n_splits))

for i in $(seq 0 $(($n_splits - 1))); do
	out_dir=../data/pdbs/splits_$i/
	mkdir "$out_dir"
	for j in $(seq 0 $(($size - 1))); do
                index=$((i * size + j))
		cp "${pdbs[index]}" "$out_dir"
	done
done

for i in $(seq $(($count - $count % $size)) $(($count - 1))); do
	echo $i
        s_idx=$((n_splits - 1))
        cp "${pdbs[i]}" "../data/pdbs/splits_$s_idx/"
done
