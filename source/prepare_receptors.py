import subprocess
import pandas as pd
import numpy as np
import argparse
import os
import json
from concurrent.futures import ProcessPoolExecutor

from utils import l2_norm
from pathlib import Path
from Bio.PDB import Structure, PDBParser
from pdbfixer import PDBFixer
from openmm.app import PDBFile
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO
from datetime import datetime

def prepare_receptor(pdb_path : Path, pocket_id, center_coords, box_sizes, out_folder='../data/docking_files') -> list[Path]:
    # Export receptor atoms
    center_x, center_y, center_z = center_coords
    size_x, size_y, size_z = box_sizes
    out_path = f'{out_folder}/{pdb_path.stem}_{pocket_id}'
    command = [
        "mk_prepare_receptor.py",
        "-i", str(pdb_path),
        "-o", out_path,
        "-p",  # Generate PDBQT file
        "-v",  # Generate Vina config
        "--box_center", str(center_x), str(center_y), str(center_z),
        "--box_size", str(size_x), str(size_y), str(size_z)
    ]
    try:
        subprocess.run(command, check=True)

    except subprocess.CalledProcessError as e:
        print(e)
        with open(f'{out_folder}/problem_pdbs.txt', 'a') as f:
            f.write(str(pdb_path) + '\n')

    return (Path(f'{out_path}.pdbqt'), Path(f'{out_path}.box.txt'))

def parse_msa_target_indices(path):
    with open(path, 'r') as f:
        lines = f.readlines()

    res = []
    for line in lines:
        left, right = line.strip().split('-')
        res.extend(list(range(int(left), int(right))))

    return res

def get_residue(structure, residue_id):
    for chain in structure.get_chains():
        for res in chain:
            if res.get_id()[1] == residue_id:
                residue = res
                break
        if residue is not None:
            break
    if residue is None:
        raise ValueError(f"Residue {residue_id} not found in structure.")
    
    return residue

def compute_pca_components(atom_coords, n_components=2):
    """Return the top PCA component vectors for the given atom coordinates."""
    coords = np.asarray(atom_coords, dtype=float)
    if coords.ndim != 2 or coords.shape[1] != 3:
        raise ValueError("atom_coords must be an array of shape (N, 3)")

    centroid = coords.mean(axis=0)
    centered = coords - centroid
    cov = np.dot(centered.T, centered) / max(centered.shape[0] - 1, 1)
    eigvals, eigvecs = np.linalg.eigh(cov)
    order = np.argsort(eigvals)[::-1]
    components = eigvecs[:, order[:n_components]].T
    return components

def compute_point_projection(point, vector, origin=None):
    point = np.asarray(point, dtype=float)
    vector = np.asarray(vector, dtype=float)
    if origin is not None:
        origin = np.asarray(origin, dtype=float)
        point = point - origin

    norm2 = np.dot(vector, vector)
    if norm2 == 0:
        raise ValueError("Projection vector must not be zero-length")

    return np.dot(point, vector) / norm2 * vector

def compute_residue_pca_projection(structure, residue_id, component=0, atom_coords=None):
    if atom_coords is None:
        atom_coords = [atom.get_coord() for atom in structure.get_atoms()]
    point = np.mean([atom.get_coord() for atom in get_residue(structure, residue_id)], axis=0)
    components = compute_pca_components(atom_coords, n_components=component + 1)
    origin = np.mean(atom_coords, axis=0)
    return compute_point_projection(point, components[component], origin=origin) + origin


def filter_pockets_pca(structure : Structure, pockets : pd.DataFrame, residues, tol=15, residue_id=5, verbose=0, include_best=False, mode='min'):
    extracellular_part_center = compute_residue_pca_projection(structure, residue_id=residue_id)
    # Mask pockets according to their distance from the centroid of selected residues 
    mask = []
    distances = []

    if verbose > 0:
        print('Pocket distance from selected residue centroid, order as in the predictions .csv:')

    if mode == 'min':
        for _, pocket in pockets.iterrows():
            min_dist = np.inf
            for residue in pocket['residue_ids'].split(' '):
                idx = int(residue[2:]) # residue ids have a pattern of A_123

                # Iterate through atoms of the residue to find the closest one
                for atom in residues[idx].get_atoms():
                    c_dist = l2_norm(atom.get_coord() - extracellular_part_center)
                    if c_dist < min_dist:
                        min_dist = c_dist
                    
            if verbose > 0:
                print(min_dist)
            distances.append(min_dist)

    else:
        distances = l2_norm(pockets[['center_x', 'center_y', 'center_z']].to_numpy(float) - extracellular_part_center) 
  
    pockets['min_surface_centroid_distance'] = distances
    mask = pockets['min_surface_centroid_distance'] < tol
    if include_best:
        # Always include the best pocket
        mask[0] = True
    if verbose > 0:
        print(mask)

    return pockets[mask]

def find_best_pockets(structure : Structure, pockets : pd.DataFrame, msa : MultipleSeqAlignment, id_to_msa_index : dict, target_indices, tol=20, mode = 'close_all',
                      verbose=1, dist_mode='min', include_best=True, k=5):
    if verbose > 0:
        print(structure.id)

    pockets = pockets.sort_values('score', ascending=False)
    if mode == 'best':
        return pockets.iloc[0]
    
    if mode == 'top_k':
        return pockets.iloc[:k]
    
    # List containing a mapping from MSA indices to last preceding sequence index
    ungapped_index = []
    idx = 0
    seq_row = id_to_msa_index[structure.id]
    for i in range(msa.get_alignment_length()):
        ungapped_index.append(idx)
        if msa[seq_row, i] != '-':
            idx += 1

    # Target indices to sequence indices
    residue_indices = []
    for i in target_indices:
        if msa[seq_row, i] != '-':
            residue_indices.append(ungapped_index[i])

    if len(residue_indices) == 0:
        # No matched residues

        # Log the protein
        if not os.path.exists(f'../temp'):
            os.makedirs('../temp', exist_ok=True)
            with open('../temp/no_match_prots.txt', 'w') as f:
                f.write(f'{structure.id}\n')

        else:
            with open(f'../temp/no_match_prots.txt', 'r') as f:
                lines = f.readlines()
                lines.append(f'{structure.id}\n')

            with open(f'../temp/no_match_prots.txt', 'w') as f:
                f.writelines(lines)

        return pockets.iloc[0]

    # Calculate a centroid from the target residues
    residues = list(structure.get_residues())
    centroid = np.average([residues[i].center_of_mass() for i in residue_indices], axis=0)

    # Mask pockets according to their distance from the centroid of selected residues 
    mask = []

    if verbose > 0:
        print('Pocket distance from selected residue centroid, order as in the predictions .csv:')

    if dist_mode == 'min':
        for _, pocket in pockets.iterrows():
            min_dist = np.inf
            for residue in pocket['residue_ids'].split(' '):
                idx = int(residue[2:]) # residue ids have a pattern of A_123
                # Iterate through atoms of the residue to find the furthest one
                try:
                    for atom in residues[idx].get_atoms():
                        c_dist = l2_norm(atom.get_coord() - centroid)
                        if c_dist < min_dist:
                            min_dist = c_dist
                except IndexError as e:
                    print(e)
                    print(structure.id)
                    print(idx)
                    print(len(list(structure.get_residues())))
                    
            #center_x, center_y, center_z = pocket['center_x'], pocket['center_y'], pocket['center_z']
            if verbose > 0:
                print(min_dist)
            mask.append(min_dist < tol)

    else:
        mask = l2_norm(pockets[['center_x', 'center_y', 'center_z']].to_numpy(float) - centroid) < tol
    
    if include_best:
        # Always include the best pocket
        mask[0] = True
    if verbose > 0:
        print(mask)
        
    # All pockets are far, return the best one
    if np.sum(mask) == 0:
        return pockets.iloc[0]

    # Return all close pockets
    if mode == 'close_all':
        return pockets[mask]
    
    # Return the pocket from the close pockets that has the highest score
    else:
        return pockets[mask].iloc[0]
    
def create_msa_index_table(msa : MultipleSeqAlignment):
    res = {}
    for i , prot in enumerate(msa):
        res[prot.id] = i

    return res

def calculate_box_size(residues, center, pocket, min_size=0, padding=2):
    max_dist = 0
    # Find the furthest atom from the center for correct box size
    for residue in pocket['residue_ids'].split(' '):
        idx = int(residue[2:]) # residue ids have a pattern of A_123
        max_dist_atom = 0
        # Iterate through atoms of the residue to find the furthest one
        for atom in residues[idx].get_atoms():
            c_dist = l2_norm(atom.get_coord() - center)
            if c_dist > max_dist_atom:
                max_dist_atom = c_dist

        if max_dist_atom > max_dist:
            max_dist = max_dist_atom
    return min(max_dist, min_size) + padding

def protonate_pdb(pdb_path : Path, ph=7, out_path=None):
    if out_path == None:
        out_path = f'../data/temp/{pdb_path.stem}_H.pdb'

    fixer = PDBFixer(str(pdb_path))
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms(seed=42)
    # subprocess.run([
    #    'reduce', '-BUILD', pdb_path, '>', out_path
    # ])
    fixer.addMissingHydrogens(ph)

    PDBFile.writeFile(fixer.topology, fixer.positions, open(out_path, 'w'))
    return Path(out_path)

def save_best_pockets(pockets, out_folder):
    os.makedirs(out_folder, exist_ok=True)
    pockets.to_csv(f'{out_folder}/filtered_pockets.csv')


def _prepare_single_pdb(task):
    pdb, args, msa, id_to_msa_index, target_indices = task
    parser = PDBParser()
    name = pdb.stem
    out_folder = os.path.join(args.out_folder, name)

    try:
        pockets = pd.read_csv(f'{args.pocket_preds_path}/{name}.pdb_predictions.csv', skipinitialspace=True)
    except FileNotFoundError as exc:
        raise FileExistsError(f'Predictions for {name} not found. Use run_p2rank.py to generate the pocket predictions first.') from exc

    h_path = protonate_pdb(pdb, args.ph)
    structure = parser.get_structure(pdb.stem, pdb)
    residues = {r.get_id()[1]: r for r in list(structure.get_residues())}
    if args.pocket_selection_mode == 'pca_close':
        best = filter_pockets_pca(structure, pockets, residues=residues, tol=args.tol, residue_id=args.pca_residue_id, include_best=args.include_best, mode=args.pocket_dist_mode)
    else:
        best = find_best_pockets(structure, pockets, msa, id_to_msa_index, target_indices, tol=args.tol,
                                mode=args.pocket_selection_mode, include_best=args.include_best,
                                k=args.k, dist_mode=args.pocket_dist_mode)

    
    save_best_pockets(best, out_folder)

    if best.shape[0] == 0:
        return []

    out = []
    if len(best.shape) > 1:
        for _, pocket in best.iterrows():
            center = pocket['center_x'], pocket['center_y'], pocket['center_z']
            center_np = np.asarray(center)
            box_size = calculate_box_size(residues, center_np, pocket, min_size=args.min_box_size, padding=args.padding)
            out.append(prepare_receptor(h_path, f'p{pocket["rank"]}', center, (box_size, box_size, box_size), out_folder=out_folder))
    else:
        center = best['center_x'], best['center_y'], best['center_z']
        center_np = np.asarray(center)
        box_size = calculate_box_size(residues, center_np, best, min_size=args.min_box_size, padding=args.padding)
        out.append(prepare_receptor(h_path, f'p{best["rank"]}', center, (box_size, box_size, box_size), out_folder=out_folder))

    return out


def prepare_receptors(args) -> list[tuple[Path, Path]]:
    pdbs = sorted(Path(args.pdbs_path).rglob('*.pdb'))
    msa = AlignIO.read(args.msa_path, format='fasta')
    id_to_msa_index = create_msa_index_table(msa)
    target_indices = parse_msa_target_indices(args.target_idxs_path)

    if not os.path.exists(args.out_folder):
        os.makedirs(args.out_folder)

    # Save config
    with open(os.path.join(args.out_folder, 'prep_receptors_config.json'), 'w') as f:
        json.dump({ k : v for k, v in vars(args).items()}, f)

    if args.n_workers and args.n_workers > 1:
        tasks = [(pdb, args, msa, id_to_msa_index, target_indices) for pdb in pdbs]
        with ProcessPoolExecutor(max_workers=args.n_workers) as executor:
            results = list(executor.map(_prepare_single_pdb, tasks))
        return [item for result in results for item in result]

    out = []
    for pdb in pdbs:
        out.extend(_prepare_single_pdb((pdb, args, msa, id_to_msa_index, target_indices)))

    return out

def add_arguments(parser : argparse.ArgumentParser):
    parser.add_argument('--pocket_preds_path', default='../data/p2rank_output', help='P2rank pocket predictions output path')
    parser.add_argument('--pdbs_path', default='../data/pdbs', help='Path to the directory where .pdb files are stored')
    parser.add_argument('--ph', default=7, type=int, help='pH for protonation')
    parser.add_argument('--msa_path', default='../data/aligned_sequences.fasta', help='MSA for the receptors, used for selecting pockets outside the membrane.')
    parser.add_argument('--target_idxs_path', default='../data/msa_index_ranges.txt', help='File containing indices to the MSA for membrane pocket selection')
    parser.add_argument('--min_box_size', default=15, type=int, help='Minimum box size for simulations')
    parser.add_argument('--padding', default=2, type=int, help='Simulation box padding')
    parser.add_argument('--tol', default=10, type=int, help='Maximum pocket center distance from the centroid determined by pocket selection')
    parser.add_argument('--pocket_selection_mode', default='pca_close', choices=['best', 'close_best', 'close_all', 'top_k', 'pca_close'], help='''
                        "best" mode simply takes the highest scoring pocket from the P2rank prediction.
                        "close_best" mode only considers the best pocket from ones that are close to the external part of the protein, determined via the indices from target_idxs_path
                        "close_all" same as above, but docks to all close pockets
                        "top_k" returns top k pockets by score
                        "pca_close" filters pockets by determining the extracellular part of the protein using the first PCA eigenvector. 
                        Center coordinates of residue with ID [pca_residue_id] are projected onto this vector, 
                        and the projected point serves as a centroid from which pocket distances are computed, and filtered according to [tol].
                        ''')
    parser.add_argument('--pca_residue_id', default=5, type=int, help='Which residue ID to consider as a reference point for the extracellular part of the protein. Relevant for pca_close mode.')
    parser.add_argument('--pocket_dist_mode', choices=['min', 'center'], default='min', help='Whether the pocket distance is calculated from the or the closest atom or the center of the pocket.')
    parser.add_argument('-k', default=5, type=int, help='If top_k mode is used, this is the k value.')
    parser.add_argument('-v', '--verbose', default=1, type=int, help='Verbosity level')
    parser.add_argument('--include_best', default=False, type=bool, help='Always include the best pocket in the selected pockets. Relevant for the "close" pocket selection mode.')
    parser.add_argument('--n_workers', default=1, type=int, help='Number of worker processes for preparing receptors in parallel. Set to 1 to disable multiprocessing.')
    parser.add_argument('-o', '--out_folder', default='../data/docking_files')
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    add_arguments(parser)
    args = parser.parse_args()
    prepare_receptors(args)
