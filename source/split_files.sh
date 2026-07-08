#!/bin/bash

set -euo pipefail

root_dir="${1%/}"
splits_dir=$2
n_splits=$3
file_pattern=$4

if [ "$#" -ne 4 ]; then
	echo "Splits files recursively into split folders, preserving the folder structure after <root_dir>"
    echo "Usage: $0 <root_dir> <splits_dir> <n_splits> <file_pattern>" >&2
	echo 'Example: ./split_tunnels.sh ../output/tunnel_results ../output/splits 64 "*.dsd"'
    exit 1
fi

mapfile -t files < <(find "$root_dir" -type f -name "$file_pattern" | sort)
file_count=${#files[@]}

echo "$file_count"

if [ "$file_count" -eq 0 ]; then
    echo "No .dsd files found." >&2
    exit 0
fi

split_size=$((file_count / n_splits))
if [ "$split_size" -eq 0 ]; then
    split_size=1
fi

echo "split size: $split_size"

current_split_size=0
split_idx=0
mkdir -p "$splits_dir"

for i in "${!files[@]}"; do
    echo "$i/$file_count"

    out_dir="$splits_dir/splits_$split_idx"
    mkdir -p "$out_dir"

    rel_path="${files[i]#"$root_dir"/}"
    target_path="$out_dir/$rel_path"
    mkdir -p "$(dirname "$target_path")"

    cp "${files[i]}" "$target_path"

    current_split_size=$((current_split_size + 1))
    if [ "$current_split_size" -eq "$split_size" ]; then
        split_idx=$((split_idx + 1))
        current_split_size=0
    fi
done
