#!/bin/bash

# Directories
PDB_DIR="../data/pdbs"
LIGAND_DIR="../data/prepared_ligands"
OUTPUT_DIR="../data/tunnel_results"
CAVER_PATH="../../caver/caver/caver.jar"
CAVER_HOME_PATH="../../caver"
CAVERDOCK_PATH="caverdock" # Assuming it's in your PATH

mkdir -p $OUTPUT_DIR

for protein in "$PDB_DIR"/*.pdb; do
    prot_name=$(basename "$protein" .pdb)
    prot_out="$OUTPUT_DIR/$prot_name"
    mkdir -p "$prot_out/caver"

    echo "--- Processing Protein: $prot_name ---"

    # 1. Run CAVER 
    # Note: You need a template config.txt where the 'pdb_file' line is updated
    #sed "s|input_pdb_file|../$protein|g" caver_template.txt > "$prot_out/caver/config.txt"
    java -Xmx4g -jar $CAVER_PATH -pdb $protein -home $CAVER_HOME_PATH -out $prot_out > caver.log

    # 2. Identify the largest tunnel (usually tunnel_1.pdb)
    TUNNEL_PDB="$prot_out/caver/out/tunnels/tunnel_1.pdb"
    break
    if [ -f "$TUNNEL_PDB" ]; then
        for ligand in "$LIGAND_DIR"/*.pdbqt; do
            lig_name=$(basename "$ligand" .pdbqt)
            dock_out="$prot_out/docking_$lig_name"
            mkdir -p "$dock_out"

            echo "Pushing $lig_name through $prot_name..."
            
            # 3. Run CaverDock
            # -i: tunnel, -l: ligand, -o: output directory
            $CAVERDOCK_PATH -i "$TUNNEL_PDB" -l "$ligand" -o "$dock_out" --full
        done
    else
        echo "Warning: No tunnel found for $prot_name"
    fi
done
