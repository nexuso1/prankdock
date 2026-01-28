import streamlit as st
import pandas as pd
import re
from pathlib import Path


# Import functions from your existing scripts
from prepare_ligands import prepare_ligands #
from prepare_receptors import prepare_receptors #
from run_p2rank import run_p2rank, create_ds_file #
from run_docking import run_docking #

# Mock Namespace to simulate argparse
class Args:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

def parse_vina_score(filepath):
    """Extracts the top docking score from a Vina output PDBQT file."""
    try:
        with open(filepath, 'r') as f:
            for line in f:
                if "REMARK VINA RESULT:" in line:
                    parts = line.split()
                    return float(parts[3])
    except Exception:
        return None
    return None

class Args:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

st.set_page_config(page_title="Vina Docking Pipeline", layout="wide")

st.title("Protein-Ligand Docking Pipeline")

tab1, tab2 = st.tabs(["Execution", "Results"])

# --- TAB 1: PIPELINE EXECUTION ---
with tab1:
    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader("General Settings")
        p2rank_path = st.text_input("P2Rank Root Directory", value="../p2rank") #
        ph_level = st.number_input("Target pH", value=7.0, step=0.1) #
        exhaustiveness = st.slider("Vina Exhaustiveness", 1, 64, 32) #
        
        st.divider()
        st.subheader("Pocket Selection")
        p_mode = st.selectbox(
            "Pocket Selection Mode",
            options=["close_all", "best", "close_best", "top_k"],
            help="'best': Highest scoring pocket. 'close_all': All pockets within tolerance."
        ) #
        k = None
        tolerance = None
        if p_mode == 'top_k':
            k = st.number_input('*k* value:', value=5)
        if p_mode == 'close_all' or 'close_best':
            tolerance = st.number_input("Distance Tolerance (Å)", value=20) #
        include_best = st.checkbox("Always include top-ranked pocket", value=True) #

    with col2:
        st.subheader("Input & Ligand Prep")
        pdb_dir = st.text_input("Directory containing .pdb files", value="../data/pdbs") #
        ligands_csv = st.text_input("Path to ligands.csv", value="../data/ligands.csv") #
        msa_path = st.text_input("MSA File (.fasta)", value="../data/aligned_sequences.fasta") #
        target_idxs = st.text_input("Target Indices File", value="../data/msa_index_ranges.txt") #
        min_box_size = st.number_input("Minimum Box Size (Å)", value=23)
        padding = st.number_input("Padding (Å)", value=2)
        st.info("Ligand Preparation Flags:")
        skip_tautomers = st.checkbox("Skip Tautomers", value=False) #
        skip_acidbase = st.checkbox("Skip Acid/Base", value=False) #

    if st.button("Start Docking Pipeline", type="primary"):
        try:
            with st.status("Running Pipeline...") as s:
                # Step 1: P2Rank
                s.update(label='Creating .ds file')
                create_ds_file(pdb_dir) #
                s.update(label='Running p2rank')
                run_p2rank(p2rank_path, '../data/pdbs.ds', '../data/p2rank_output') #
                
                # Step 2: Prepare Receptors
                s.update(label='Preparing receptors')
                rec_args = Args(pdbs_path=pdb_dir, msa_path=msa_path, target_idxs_path=target_idxs, 
                                ph=ph_level, pocket_preds_path='../data/p2rank_output', 
                                tol=tolerance, pocket_selection_mode=p_mode, k=k,
                                include_best=include_best, verbose=0,
                                min_box_size=min_box_size, padding=padding) #

                prepare_receptors(rec_args) #

                # Step 3: Prepare Ligands
                s.update(label='Preparing ligands')
                lig_args = Args(ligands_path=ligands_csv, ph=ph_level, 
                                skip_tautomers=skip_tautomers, skip_acidbase=skip_acidbase) #
                prepare_ligands(lig_args) #

                # Step 4: Docking
                s.update(label="Docking")
                dock_args = Args(dock_files_path='../data/docking_files', 
                                 ligands_path='../data/prepared_ligands', exhaustiveness=exhaustiveness) #
                
                run_docking(dock_args) #
                s.update(label="Complete!", state="complete")
            st.success("All steps finished.")
        except Exception as e:
            st.error(f"Error: {e}")

# --- (Inside Tab 2: Results Gallery) ---
with tab2:
    output_base = Path("../output")
    
    if not output_base.exists():
        st.info("No results found. Run the pipeline first.")
    else:
        receptor_dirs = [d for d in output_base.iterdir() if d.is_dir()]
        
        if not receptor_dirs:
            st.warning("Output directory is empty.")
        else:
            sel_rec_dir = st.selectbox("Select Receptor", receptor_dirs, format_func=lambda x: x.name)
            
            # --- POCKET PARSING LOGIC ---
            # Result folders are named like 'protein_p1'. We extract 'protein' and '1'.
            folder_name = sel_rec_dir.name
            match = re.search(r'^(.*)_H_p(\d+)$', folder_name)
            
            if match:
                pdb_stem = match.group(1)
                p_id = match.group(2)
                # Fetch residue IDs for this pocket from P2Rank output

            # ----------------------------

            receptor_pdbqt = Path(f"../data/docking_files/{folder_name}.pdbqt")
            result_files = list(sel_rec_dir.glob("*_out.pdbqt"))
            data = [{"Ligand": f.name.replace("_out.pdbqt", ""), "Affinity": parse_vina_score(f), "Path": f} for f in result_files]
            
            if data:
                df = pd.DataFrame(data).sort_values(by="Affinity")
                st.dataframe(df[["Ligand", "Affinity"]], width='stretch', hide_index=True)