#!/bin/bash
count=0
declare qts=()
declare boxes=()
declare configs=()
for file in $(ls ./docking_files/*.pdbqt); do
        qts+=($file)
        boxes+=(${file%.pdbqt}.box.pdb)
        configs+=(${file%.pdbqt}.box.txt)
	count=$((count+1))
done

echo $count
echo ${qts[0]}
echo ${boxes[0]}
echo ${configs[0]}

# Modify this based on desired number of jobs
n_splits=20
size=$((count / $n_splits))

for i in $(seq 0 $(($n_splits - 1))); do
	mkdir ./splits_$i
	for j in $(seq 0 $size); do
                index=$((i * size + j))
		fname=
		cp "${qts[index]}" "./splits_$i/"
                cp "${boxes[index]}" "./splits_$i/"
		cp "${configs[index]}" "./splits_$i/"
	done
done

for i in $(seq $(($count - $count % $size)) $(($count - 1))); do
	echo $i
        s_idx=$((n_splits - 1))
        cp "${qts[i]}" "./splits_$s_idx/"
        cp "${boxes[i]}" "./splits_$s_idx/"
	cp "${configs[i]}" "./splits_$s_idx/"
done
