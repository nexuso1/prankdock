# PrankDock - Automated pocket prediction and docking pipeline



## Setup

### 1. Environmnent setup
- Download and install the Anaconda package manager, and create a new virtual environment for prankdock using the `.yaml` file.    
```
conda env create -f environment.yaml
```
- Activate the environment using
```
conda activate dock
```
- Run the setup script
```
cd source
./setup.sh
```

### 2. PDB Downloads
If you wish to download PDBs via `download_pdbs.py`, you need to obtain an API key for AlphaFold Database. Create a `.txt` file containing the key at this path: `data/key.txt`.

## Usage
- NOTE: The scripts mentioned in this section expect to be run from the `source/` directory, so first run
```
cd source
```

- If you're unsure about anything, you can run any script with the `-h` or `--help` option to see argument explanations.
1. Download PDBs from AFDB via `source/download_pdbs.py`. 
    - This script requires a `.csv`containing the Uniprot IDs of requested receptors, inside a column named `uniprot_id`.
    - Default path is `data/prots.csv`
    - If you already have PDBs downloaded, place them into `data/pdbs`.
2. Predict pockets with P2rank using `source/run_p2rank.py`.
    -  By default, the script gathers the paths of all `.pdb` files in `data/pdbs`, runs P2rank, and deposits the results into `data/p2rank_output`
4. Prepare ligands using `source/prepare_ligands.py`.
    - This script requires a `.csv` file containing the following fields: `name` and `smiles`. `name` is the name of the ligand, and `smiles` is its smiles string.
    - Default location of this file is `data/ligands.csv`
    - If you wish to skip isomers, use the `--skip_tautomer` and `--skip_acidbase` arguments

5. Prepare receptors using `source/prepare_receptors.py`

    - This script handles receptor preparation, including addition of hydrogen atoms and docking pocket selection.
        - Pocket selection for docking has a few modes, the default one is `pca_close`. It works as follows:
            - First, we compute the first PCA component for the given protein, giving us a good central axis estimate
            - Then, a centroid is computed from the projection of the atomic coordinates of a residue with a specifc sequential index (residue 5 by default) onto the central axis.
                - Since the N-terminal should be near the external part of the membrane, we choose a residue that is close to the beginning of the sequence
                - The first residue can potentially move a lot, so we choose a residue with a slightly higher sequential index
            - For each pocket, it computes the distance between the closest atom in the pocket and the centroid
            - If the pocket is close enough (<= 10 A), it is considered for docking
        - use `python run_docking.py -h` to get information about the other modes
    - Pocket size is computed automatically according to the p2rank pocket prediction
    - Each protonated receptor-pocket pair is then passed to `mk_prepare_receptor.py`
    - Prepared receptors and vina configs are stored in `data/docking_files/<protein_id>/po` by default
3. Dock ligands into predicted pockets with `source/run_docking.py`. 

    - After the receptors are prepared, it runs the AutoDock Vina 4 using the Vina forcefield on each combination of receptor and ligand.
    - Outputs are deposited in the `output/` folder
        - each pocket has its own subfolder (e.g. `output/output/A0A067XG43_p1` for pocket with rank 1 from the P2rank prediction) containing the `.pdbqt` files of docked ligands

## Example
- assuming you already have downloaded the PDBs

```
conda activate dock
cd source
python run_p2rank.py
python prepare_ligands.py
python prepare_receptors.py
python run_docking.py
```

## Gui
- You can run the graphical interface locally from the source folder using
```
streamlit run gui.py
```
Your browser should open a new tab where you can interact with the docking pipeline and view results.
<img width="1821" height="834" alt="image" src="https://github.com/user-attachments/assets/24e1b18e-969a-42db-a44c-ba8e349f9795" />
<img width="1789" height="569" alt="image" src="https://github.com/user-attachments/assets/449ea279-2596-4e8d-9c72-ba4ccc2c64e2" />


# Tunnel Pipeline
For each protein, computes possible tunnels and pushes all ligands through them. 
Currently assumes you are running on a machine that is able to use MPI to speed

## Requirements
- Apptainer
    - https://apptainer.org/
- Conda
    - https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html

## Usage
1. Run the CaverDock container using Apptainer/Singularity
```
apptainer shell caverdock-1.2.sif
```
2. Run the `tunnel_pipeline.sh` script
```sh
# Example usage
# Assumes you already ran p2rank and prepared both ligands and receptors.

cd source

# pockets 1-3 = dock into pockets ranked 1 through 3 (1 = highest score)
# tunnels 1-3 = use tunnels 1 through 3 (1 = highest score)
# delta 0.3 = Analyze ligand energy every 0.3 Ångström while going through the tunnel
# exhaustiveness = Vina exhaustiveness (higher values MASSIVELY increase computation time)
./tunnel_pipeline.sh --pockets 1-3 --tunnels 1-3 --exhaustiveness 4 --delta 0.3 --o ../output/tunnel_results --caver_jar ../caver/caver.jar
```