#!/bin/bash

# --- Argument Parsing ---
PDB_DIR="../data/pdbs"
MPI_PROCS=1
POCKET_START=1
POCKET_END=1
TUNNEL_START=1
TUNNEL_END=1
DELTA=0.5
EXHAUSTIVENESS=16
SHELL_DEPTH=2.5
CLUSTERING_THRESHOLD=4.5
OUTPUT_DIR="../output/tunnel_results"
CAVER_PATH="../../caver/caver.jar"
while [[ $# -gt 0 ]]; do
    case "$1" in
        -p|--pockets)
            IFS='-' read -r POCKET_START POCKET_END <<< "$2"
            shift 2 ;;
        -t|--tunnels)
            IFS='-' read -r TUNNEL_START TUNNEL_END <<< "$2"
            shift 2 ;;
        --delta)
            DELTA="$2"
            shift 2 ;;
        -e|--exhaustiveness)
            EXHAUSTIVENESS="$2"
            shift 2 ;;
        --mpi_procs)
            MPI_PROCS="$2"
            shift 2 ;;

        -o|--output_dir)
            OUTPUT_DIR="$2"
            shift 2 ;;
        -i|--pdb_dir)
            PDB_DIR="$2"
            shift 2 ;;
        -d|--shell_depth)
            SHELL_DEPTH="$2"
            shift 2 ;;
        -ct|--clustering_threshold)
            CLUSTERING_THRESHOLD="$2"
            shift 2 ;;
        --caver_jar)
            CAVER_PATH="$2"
            shift 2 ;;
        *)
            if [ -z "$PDB_DIR" ]; then
                PDB_DIR="$1"
            elif [ -z "$MPI_PROCS" ]; then
                MPI_PROCS="$1"
            fi
            shift ;;
    esac
done

echo "Running on data in $PDB_DIR with $MPI_PROCS MPI processes"
echo "Pockets: $POCKET_START-$POCKET_END | Tunnels: $TUNNEL_START-$TUNNEL_END"

cat <<EOF > $OUTPUT_DIR/config.txt
argument value
pockets $POCKET_START-$POCKET_END
tunnels $TUNNEL_START-$TUNNEL_END
delta $DELTA
exhaustiveness $EXHAUSTIVENESS
shell_depth $SHELL_DEPTH
clustering_threshold $CLUSTERING_THRESHOLD
EOF

# Configuration
PDB_DIR="../data/pdbs"
PREP_RECEPTOR_DIR="../data/docking_files"
LIGAND_DIR="../data/prepared_ligands"
CSV_DIR="../data/p2rank_output"

CAVER_HOME_PATH=$(dirname "$CAVER_PATH")

mkdir -p $OUTPUT_DIR

for pdb_path in "$PDB_DIR"/*.pdb; do
    prot_name=$(basename "$pdb_path" .pdb)
    echo "--- Processing $prot_name ---"
    csv_file="$CSV_DIR/${prot_name}.pdb_predictions.csv"

    if [ ! -f "$csv_file" ]; then
        echo "Skipping $prot_name: No CSV found at $csv_file"
        continue
    fi

    # --- Pocket loop ---
    for pocket_idx in $(seq "$POCKET_START" "$POCKET_END"); do
        pocket_name="pocket${pocket_idx}"
        echo "=== Pocket $pocket_idx ($pocket_name) for $prot_name ==="

        coords=$(awk -F',' -v pocket="$pocket_name" '$1 ~ pocket {print $7, $8, $9}' "$csv_file")
        read -r cx cy cz <<< "$coords"
        echo "Starting coords for $prot_name ($pocket_name): $coords"

        if [ -z "$cx" ]; then
            echo "Warning: Could not find $pocket_name coordinates for $prot_name; skipping this pocket"
            continue
        fi

        # Each pocket gets its own output subdirectory
        prot_out_dir="$OUTPUT_DIR/$prot_name/pocket${pocket_idx}"
        mkdir -p "$prot_out_dir"
        config_path="$prot_out_dir/config.txt"
        cat <<EOF > "$config_path"
starting_point_coordinates $cx $cy $cz
frame_clustering_threshold $CLUSTERING_THRESHOLD
shell_depth $SHELL_DEPTH
save_dynamics_visualization yes
EOF
        # Copy PDB once per pocket directory (CAVER reads from the -pdb dir)
        cp "$pdb_path" "$prot_out_dir/$prot_name.pdb"

        run_caver() {
            local conf="$1"
            rm -rf "$prot_out_dir/data"
            echo "Running CAVER for $prot_name (pocket${pocket_idx}) with config $(basename "$conf"), saving to $prot_out_dir"
            java -Xmx4g -jar "$CAVER_PATH" -conf "$conf" -pdb "$prot_out_dir" -home "$CAVER_HOME_PATH" -cp "$CAVER_HOME_PATH"/lib -out "$prot_out_dir" > "$prot_out_dir/caver.log" 2>&1
        }

        run_discretizer() {
            echo "Discretizing tunnel for $prot_name (pocket${pocket_idx}, tunnel${tunnel_idx})"
            discretizer -f "$TUNNEL" -o "$discr_tunnel" --delta "$DELTA" > "$prot_out_dir/discretizer_tunnel${tunnel_idx}.log" 2>&1
        }

        run_caver "$config_path"

        # --- Tunnel loop ---
        for tunnel_idx in $(seq "$TUNNEL_START" "$TUNNEL_END"); do
            # CAVER zero-pads cluster IDs to 3 digits: 001, 002, ...
            tunnel_id=$(printf "%03d" "$tunnel_idx")
            TUNNEL="$prot_out_dir/data/clusters_timeless/tun_cl_${tunnel_id}_1.pdb"

            if [ ! -f "$TUNNEL" ]; then
                echo "Tunnel $tunnel_idx not found for $prot_name (pocket${pocket_idx}); skipping"
                continue
            fi

            echo "--- Tunnel $tunnel_idx for $prot_name (pocket${pocket_idx}) ---"

            discr_tunnel="$prot_out_dir/data/${tunnel_idx}.dsd"
            prep_receptor="$PREP_RECEPTOR_DIR"/"$prot_name"_H_p1.pdbqt

            if ! run_discretizer; then
                if grep -q 'AssertionError' "$prot_out_dir/discretizer_tunnel${tunnel_idx}.log"; then
                    echo "AssertionError detected. Rerunning CAVER with awvd no."
                    alt_config_path="$prot_out_dir/config-awvd-no.txt"
                    cat <<EOF > "$alt_config_path"
starting_point_coordinates $cx $cy $cz
frame_clustering_threshold $CLUSTERING_THRESHOLD
shell_depth $SHELL_DEPTH
save_dynamics_visualization yes
awvd no
EOF
                    run_caver "$alt_config_path"

                    if [ ! -f "$TUNNEL" ]; then
                        echo "No tunnel $tunnel_idx found after CAVER rerun with awvd no; skipping"
                        continue
                    fi

                    echo "Retrying discretizer after CAVER rerun"
                    if ! run_discretizer; then
                        echo "Discretizer still failed; see $prot_out_dir/discretizer_tunnel${tunnel_idx}.log"
                        continue
                    fi
                else
                    echo "Discretizer failed; see $prot_out_dir/discretizer_tunnel${tunnel_idx}.log"
                    continue
                fi
            fi

            # --- Ligand docking loop ---
            for ligand in "$LIGAND_DIR"/*.pdbqt; do
                l_name=$(basename "$ligand" .pdbqt)
                dock_out_dir="$prot_out_dir/docking_${l_name}_tunnel${tunnel_idx}"
                mkdir -p "$dock_out_dir"

                result_file="$dock_out_dir/${l_name}_tunnel${tunnel_idx}-ub.pdbqt"
                if [ -f "$result_file" ]; then
                    echo "Skipping already completed: $l_name into $prot_name pocket${pocket_idx} tunnel${tunnel_idx}"
                    continue
                fi

                dock_conf="$dock_out_dir/cd_tunnel${tunnel_idx}.conf"
                cd-prepareconf -r "$prep_receptor" -l "$ligand" -t "$discr_tunnel" --seed 42 --exhaustiveness "$EXHAUSTIVENESS" > "$dock_conf"
                echo "Docking $l_name into $prot_name (pocket${pocket_idx}) tunnel $tunnel_idx ($discr_tunnel)"

                if [ "$MPI_PROCS" -eq 1 ]; then
                    caverdock --config "$dock_conf" --out "$dock_out_dir/${l_name}_tunnel${tunnel_idx}" > "$dock_out_dir/log.txt"
                else
                    mpirun.openmpi -v -np "$MPI_PROCS" caverdock --config "$dock_conf" --out "$dock_out_dir/${l_name}_tunnel${tunnel_idx}" > "$dock_out_dir/log.txt"
                fi
                energy_dat="$dock_out_dir/tunnel${tunnel_idx}_energy_profile.dat"
                energy_pdf="$dock_out_dir/tunnel${tunnel_idx}_energy_profile.png"
                cd-energyprofile -d "$discr_tunnel" -t "$result_file" -s 0 > "$energy_dat"
                gnuplot -e "set terminal png; set output '$energy_pdf'; set xlabel 'distance'; set ylabel 'energy'; plot '$energy_dat' u 1:4 w l t 'upper-bound', '' u 1:6 w l t 'lower-bound'"
            done

        done  # tunnel loop
    done  # pocket loop
done  # protein loop