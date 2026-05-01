#!/bin/bash

# Configuration
PDB_DIR="../data/pdbs"
LIGAND_DIR="../data/prepared_ligands"
CSV_DIR="../data/p2rank_output"  # Folder containing the CSVs
OUTPUT_DIR="../data/tunnel_results"
CAVER_PATH="../../caver/caver.jar"

mkdir -p $OUTPUT_DIR

for pdb_path in "$PDB_DIR"/*.pdb; do
    prot_name=$(basename "$pdb_path" .pdb)
    csv_file="$CSV_DIR/${prot_name}.pdb_predictions.csv"
    
    if [ ! -f "$csv_file" ]; then
        echo "Skipping $prot_name: No CSV found at $csv_file"
        continue
    fi

    echo "--- Processing $prot_name ---"

    # Extract coordinates for pocket1 from CSV
    # Assumes headers: name, rank, score... center_x(7), center_y(8), center_z(9)
    # We use awk to find the line where the first column is 'pocket1'
    coords=$(awk -F',' '$1 ~ /pocket1/ {print $7, $8, $9}' "$csv_file")
    read -r cx cy cz <<< "$coords"
    echo "Starting coords for $prot_name: $coords"
    if [ -z "$cx" ]; then
        echo "Error: Could not find pocket1 coordinates for $prot_name"
        continue
    fi

    # Create protein-specific CAVER config
    config_path="$OUTPUT_DIR/$prot_name/config.txt"
    cat <<EOF > "$config_path"
input_pdb_file $pdb_path
output_directory $OUTPUT_DIR
starting_point_coordinates $cx $cy $cz
shell_radius 3.0
shell_depth 4.0
frame_clustering_threshold 3.5
EOF

    # Run CAVER
    java -Xmx4g -jar "$CAVER_PATH" config.txt > $OUTPUT_DIR/$prot_name/caver.log


    # Run CaverDock for the best tunnel (tunnel_1.pdb)
    TUNNEL="$OUTPUT_DIR/$prot_name/caver/out/tunnels/tunnel_1.pdb"
    
    if [ -f "$TUNNEL" ]; then
        for ligand in "$LIGAND_DIR"/*.pdbqt; do
            l_name=$(basename "$ligand" .pdbqt)
            dock_out="$OUTPUT_DIR/$prot_name/docking_$l_name"
            mkdir -p "$dock_out"
            
            echo "Docking $l_name into $prot_name tunnel..."
            caverdock -i "$TUNNEL" -l "$ligand" -o "$dock_out" --full > "$dock_out/log.txt"
        done
    fi
done