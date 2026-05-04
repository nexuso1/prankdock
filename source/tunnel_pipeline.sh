#!/bin/bash

SPLIT=$1
MPI_PROCS=$2
echo "Running split $SPLIT with $MPI_PROCS MPI processes"

# Configuration
PDB_DIR="../data/pdbs"
PREP_RECEPTOR_DIR="../data/docking_files"
LIGAND_DIR="../data/prepared_ligands"
CSV_DIR="../data/p2rank_output"  # Folder containing the CSVs
OUTPUT_DIR="../data/tunnel_results"
CAVER_PATH="../../caver/caver.jar"
CAVER_HOME_PATH="../../caver"
mkdir -p $OUTPUT_DIR

for pdb_path in "$PDB_DIR"/splits_"$SPLIT"/*.pdb; do
    prot_name=$(basename "$pdb_path" .pdb)
    echo "--- Processing $prot_name ---"
    csv_file="$CSV_DIR/${prot_name}.pdb_predictions.csv"

    if [ ! -f "$csv_file" ]; then
        echo "Skipping $prot_name: No CSV found at $csv_file"
        continue
    fi

    

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
    prot_out_dir="$OUTPUT_DIR/$prot_name"
    mkdir -p "$prot_out_dir"
    config_path="$prot_out_dir/config.txt"
    cat <<EOF > "$config_path"
starting_point_coordinates $cx $cy $cz
frame_clustering_threshold 4.5
shell_depth 2.5
EOF
    cp "$pdb_path" "$prot_out_dir/$prot_name.pdb"

    run_caver() {
        local conf="$1"
        rm -rf "$prot_out_dir/data"
        echo "Running CAVER for $prot_name with config $(basename "$conf")"
        java -Xmx4g -jar "$CAVER_PATH" -conf "$conf" -pdb "$prot_out_dir" -home "$CAVER_HOME_PATH" -cp "$CAVER_HOME_PATH"/lib -out "$prot_out_dir" > "$prot_out_dir/caver.log" 2>&1
    }

    run_discretizer() {
        echo "Discretizing tunnel for $prot_name"
        discretizer -f "$TUNNEL" -o "$discr_tunnel" --delta 0.5 > "$prot_out_dir/discretizer.log" 2>&1
    }

    run_caver "$config_path"

    # Run CaverDock for the best tunnel (1.pdb)
    TUNNEL="$prot_out_dir/data/clusters_timeless/tun_cl_001_1.pdb"
    
    if [ -f "$TUNNEL" ]; then
        discr_tunnel="$prot_out_dir/data/1.dsd"
        prep_receptor="$PREP_RECEPTOR_DIR"/"$prot_name"_H_p1.pdbqt
        
        if ! run_discretizer; then
            if grep -q 'AssertionError' "$prot_out_dir/discretizer.log"; then
                echo "AssertionError detected during discretizer. Rerunning CAVER with awvd no."
                alt_config_path="$prot_out_dir/config-awvd-no.txt"
                cat <<EOF > "$alt_config_path"
starting_point_coordinates $cx $cy $cz
frame_clustering_threshold 4.5
shell_depth 2.5
awvd no
EOF
                run_caver "$alt_config_path"
                if [ ! -f "$TUNNEL" ]; then
                    echo "No tunnel found after rerunning CAVER with awvd no; skipping $prot_name"
                    continue
                fi
                echo "Retrying discretizer after CAVER rerun"
                if ! run_discretizer; then
                    echo "Discretizer still failed after awvd no retry; see $prot_out_dir/discretizer.log"
                    continue
                fi
            else
                echo "Discretizer failed for $prot_name; see $prot_out_dir/discretizer.log"
                continue
            fi
        fi

        # Dock all ligands
        for ligand in "$LIGAND_DIR"/*.pdbqt; do
            l_name=$(basename "$ligand" .pdbqt)
            dock_out_dir="$prot_out_dir/docking_$l_name"
            mkdir -p "$dock_out_dir"

            # Skip already completed ligand dockings
            if [ -f "$dock_out_dir"/"$l_name"_tunnel1-ub.pdbqt ]; then
                continue
            fi

            dock_conf="$dock_out_dir"/cd_tunnel1.conf
            cd-prepareconf -r "$prep_receptor" -l "$ligand" -t "$discr_tunnel" --seed 42 --exhaustiveness 16 > "$dock_conf"
            echo "Docking $l_name into $prot_name ($prep_receptor) tunnel 1 ($discr_tunnel)"

            # In order to run deterministically, we need to use -np 2 and a seed. But this takes forever.
            # mpirun.openmpi is used instead of mpirun because mpirun didnt work inside the caverdock container while inside
            # a job on the elmo cluster
            mpirun.openmpi -v -np "$MPI_PROCS" caverdock --config "$dock_conf" --out "$dock_out_dir"/"$l_name"_tunnel1 > "$dock_out_dir/log.txt"
            cd-energyprofile -d "$discr_tunnel" -t "$dock_out_dir"/"$l_name"_tunnel1-ub.pdbqt -s 0 > "$dock_out_dir"/tunnel1_energy_profile.dat
        done
    fi
done
