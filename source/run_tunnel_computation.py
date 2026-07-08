#!/usr/bin/env python3
import argparse
import csv
import re
import shutil
import subprocess
import pandas as pd
import numpy as np
import glob
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Sequence
from prepare_receptors import compute_residue_pca_projection
from Bio.PDB import PDBParser
from utils import l2_norm

@dataclass
class TunnelingConfig:
    pdb_dir: str | Path = "../data/pdbs"
    pocket_predictions_dir: str | Path = "../data/docking_files/filtered_receptors"
    output_dir: str | Path = "../output/tunnel_results"
    caver_jar: str | Path = "../caver_3.0/caver/caver.jar"
    prep_receptor_dir: str | Path = "../data/docking_files"
    cd_image_path:str |Path = "../caverdock-1.2.sif"
    pockets: str = "all"
    tunnels: str = "all"
    delta: float = 0.5
    shell_depth: float = 2.5
    clustering_threshold: float = 4.5
    centroid_res_id: int = 5
    tol : float = 5.0
    filter_tunnels: bool = True
    compute_tunnel_residues: bool = True

    def __str__(self):
        buf = ['argument,value']
        for attr, val in vars(self).items():
            buf.append(f'{attr},{str(val)}')

        return '\n'.join(buf)

def create_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run CAVER tunnel computation for proteins and pockets")
    parser.add_argument("--pdb_dir", default="../data/protonated_pdbs", help="Folder containing input PDB files")
    parser.add_argument("--preds_dir", default="../data/docking_files/filtered_receptors", help="Folder containing p2rank pocket coordinate CSV files")
    parser.add_argument("--output_dir", default="../output/tunnel_results", help="Where tunnel outputs will be written")
    parser.add_argument("--caver_jar", default="../caver_3.0/caver/caver.jar", help="Path to the CAVER jar file")
    parser.add_argument("--prep_receptor_dir", default="../data/docking_files", help="Prepared receptor directory (kept for compatibility)")
    parser.add_argument("--pockets", default="all", help="Pockets to process. Either 'all' or a range (for example 1-3)")
    parser.add_argument("--tunnels", default="all", help="Tunnel to process. Either 'all' or a range (for example 1-3)")
    parser.add_argument("--delta", type=float, default=0.3, help="Discretizer delta. 0.3A is the Caverdock default.")
    parser.add_argument("--shell_depth", type=float, default=4, help="CAVER shell depth")
    parser.add_argument("--clustering_threshold", type=float, default=5, help="CAVER clustering threshold")
    parser.add_argument("--compute_tunnel_residues", type=bool, default=False, help="Compute tunnel residues")
    parser.add_argument("--tol", type=float, default=10, help="Tunnel filtering surface distance tolerance, in Angstroms.")
    parser.add_argument("--centroid_res_id", type=int, default=5, help="Index of the resiude from which a surface centroid will be computed.")
    parser.add_argument("--cd_image", type=str, default="../caverdock-1.2.sif", help="Path to CaverDock Apptainer image")
    return parser

def parse_range(value: str) -> tuple[int, int] | None:
    normalized = value.strip()
    if normalized.lower() == "all":
        return None

    parts = normalized.split("-", 1)
    start = int(parts[0])
    end = int(parts[1]) if len(parts) > 1 else start
    if end < start:
        raise ValueError("end must be greater than or equal to start")
    return start, end


def discover_pocket_indices(csv_path: Path) -> list[int]:
    indices: list[int] = []
    with csv_path.open(newline="") as handle:
        reader = csv.reader(handle, skipinitialspace=True)
        next(reader, None)
        for row in reader:
            if not row:
                continue
            match = re.match(r"pocket(\d+)", row[0].strip())
            if match:
                indices.append(int(match.group(1)))
    return sorted(set(indices))

def resolve_path(path_value: str | Path) -> Path:
    path = Path(path_value)
    return path.resolve()

def read_pocket_coordinates(csv_path: Path, pocket_name: str) -> tuple[float, float, float] | None:
    with csv_path.open(newline="") as handle:
        reader = csv.reader(handle, skipinitialspace=True)
        next(reader, None)
        for row in reader:
            if not row:
                continue
            if row[0].strip() == pocket_name:
                try:
                    return float(row[6]), float(row[7]), float(row[8])
                except (IndexError, ValueError):
                    return None
    return None

def get_tunnel_endpoints(pdb_filename):
    first_coords = None
    last_coords = None

    with open(pdb_filename, 'r') as file:
        for line in file:
            # Only process lines starting with ATOM or HETATM
            if line.startswith("ATOM") or line.startswith("HETATM"):
                # Split the line by any whitespace
                parts = line.split()

                # In a standard ATOM line split, indices are:
                # parts[6] -> X, parts[7] -> Y, parts[8] -> Z
                x = float(parts[6])
                y = float(parts[7])
                z = float(parts[8])
                
                current_coords = (x, y, z)
                
                if first_coords is None:
                    first_coords = current_coords
                
                last_coords = current_coords

    return first_coords, last_coords

def write_config(path: Path, center_coords : tuple[float], clustering_threshold: float, shell_depth: float, awvd: bool = False, compute_tunnel_residues=False) -> None:
    cx, cy, cz = center_coords
    lines = [
        f"starting_point_coordinates {cx} {cy} {cz}",
        f"frame_clustering_threshold {clustering_threshold}",
        f"shell_depth {shell_depth}",
        f"compute_tunnel_residues {'yes' if compute_tunnel_residues else 'no'}",
        "residue_contact_distance 2.0",
        "save_dynamics_visualization yes",
        "seed 42",
        "swap no",
        "awvd no"
        "automatic shell radius yes"
    ]
    path.write_text("\n".join(lines) + "\n")


def run_command(command: Sequence[str], log_path: Path, cwd: Optional[Path] = None) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w") as handle:
        completed = subprocess.run(list(command), cwd=cwd, stdout=handle, stderr=subprocess.STDOUT, check=False)
    if completed.returncode != 0:
        raise subprocess.CalledProcessError(completed.returncode, command)


def run_caver(conf_path: Path, prot_out_dir: Path, pdb_path: Path, caver_jar: Path) -> None:
    caver_home = caver_jar.parent
    command = [
        "java",
        "-Xmx4g",
        "-jar",
        str(caver_jar),
        "-conf",
        str(conf_path),
        "-pdb",
        str(prot_out_dir),
        "-home",
        str(caver_home),
        "-cp",
        str(caver_home / "lib"),
        "-out",
        str(prot_out_dir),
    ]
    shutil.copy2(pdb_path, prot_out_dir / pdb_path.name)
    run_command(command, prot_out_dir / "caver.log")

def save_tunnel_config(output_dir: Path, config) -> None:
    config_path = output_dir / "tunneling_config.csv"
    config_path.write_text(str(config))

def get_pocket_indices(range: str, pocket_preds : pd.DataFrame) -> list[int]:
    if range.lower() == "all":
        return pocket_preds.index

    pocket_range = parse_range(range)
    if pocket_range is None:
        raise ValueError("Unable to parse pocket range")
    indices = set(list(range(pocket_range[0], pocket_range[1] + 1)))
    mask = pocket_preds.name.apply(lambda x: x.strip().endswith(idx) for idx in indices)
    return pocket_preds.index[mask]

def get_tunnel_indices(tunnels_arg: str, prot_out_dir: Path) -> list[int]:
    if tunnels_arg.lower() == "all":
        tunnel_dir = prot_out_dir / "data" / "clusters_timeless"
        if not tunnel_dir.exists():
            return []
        available_tunnels = sorted(tunnel_dir.glob("tun_cl_*_1.pdb"))
        return sorted({int(path.stem.split("_")[2]) for path in available_tunnels})

    tunnel_range = parse_range(tunnels_arg)
    if tunnel_range is None:
        raise ValueError("Unable to parse tunnel range")
    return list(range(tunnel_range[0], tunnel_range[1] + 1))

def discretize_tunnel(
    tunnel_path: Path,
    output_path: Path,
    log_path: Path,
    delta: float,
    cd_image_path: Path,
) -> Path | None:

    if not tunnel_path.exists():
        print(f"Tunnel not found; skipping")
        return None

    discretizer_cmd = ["apptainer", "exec",  str(cd_image_path), "discretizer", "--file", str(tunnel_path), "-o", str(output_path), "--delta", str(delta)]

    try:
        run_command(discretizer_cmd, log_path)
    except subprocess.CalledProcessError:
        print(f"Discretizer failed; see {log_path}")
        return None

    print(f"Saved discretized tunnel to {output_path}")
    return output_path

def filter_tunnels(tunnel_endpoints, pdb_path, centroid_residue_id, tol=5):
    parser = PDBParser()
    structure = parser.get_structure(pdb_path.stem, pdb_path)
    extracellular_part_center = compute_residue_pca_projection(structure, residue_id=centroid_residue_id)
    res = []
    for tunnel, endpoints in tunnel_endpoints.items():
        min_dist = np.inf
        for point in endpoints:

            # Iterate through atoms of the residue to find the closest one
                c_dist = l2_norm(point - extracellular_part_center)
                if c_dist < min_dist:
                    min_dist = c_dist

        if min_dist < tol:
            res.append(tunnel)

    return res
 
def process_pocket(
    pocket_idx: int,
    prot_name: str,
    pdb_path: Path,
    pocket_center_coords: tuple[float],
    config: TunnelingConfig,
) -> list[Path]:
    pocket_name = f"pocket{pocket_idx}"
    if pocket_center_coords is None:
        print(f"Warning: could not find {pocket_name} coordinates for {prot_name}; skipping")
        return []

    cx, cy, cz = pocket_center_coords
    prot_out_dir = Path(config.output_dir) / prot_name / pocket_name
    prot_out_dir.mkdir(parents=True, exist_ok=True)
    config_file = prot_out_dir / "config.txt"
    write_config(config_file, pocket_center_coords, config.clustering_threshold, config.shell_depth, compute_tunnel_residues=config.compute_tunnel_residues)
    shutil.copy2(pdb_path, prot_out_dir)
    shutil.rmtree(prot_out_dir / "data", ignore_errors=True)

    print(f"=== Pocket {pocket_idx} ({pocket_name}) for {prot_name} ===")
    print(f"Starting coords for {prot_name} ({pocket_name}): {cx}, {cy}, {cz}")

    try:
        run_caver(config_file, prot_out_dir, pdb_path, Path(config.caver_jar))
    except subprocess.CalledProcessError as exc:
        print(f"CAVER exited with code {exc.returncode}; see {prot_out_dir / 'caver.log'}")
        return []

    created_tunnels: list[Path] = []
    tunnel_endpoints = { Path(path) : get_tunnel_endpoints(path) for path in glob.glob(str(prot_out_dir / 'data' / 'clusters' / '*.pdb')) }
    print(f"Number of tunnels: {len(tunnel_endpoints.keys())}")
    for tunnel_path in filter_tunnels(tunnel_endpoints, pdb_path, config.centroid_res_id, config.tol):
        tunnel_id = tunnel_path.stem.removesuffix('.pdb').split('_')[-1]
        discr_tunnel = prot_out_dir / "data" / f"tunnel_{tunnel_id}.dsd"
        log_path = prot_out_dir / f"discretizer_tunnel{tunnel_id}.log"
        print(f"--- Tunnel {tunnel_id} for {prot_name} ({pocket_name}) ---")
        result_path = discretize_tunnel(
            tunnel_path=tunnel_path,
            output_path=discr_tunnel,
            log_path=log_path,
            delta=config.delta,
            cd_image_path=config.cd_image_path,
        )
        if result_path is not None:
            created_tunnels.append(result_path)
    print(f"Number of tunnels post-filtering: {len(created_tunnels)}")
    return created_tunnels

def get_predictions_csv(prot_name, preds_dir):
    preds_dir = Path(preds_dir)
    csv_file = preds_dir / f"{prot_name}.pdb_predictions.csv"
    if not csv_file.exists():
        csv_file = preds_dir / prot_name / 'filtered_pockets.csv'
    
    if not csv_file.exists():
        print(f"Skipping {prot_name}: no p2rank CSV found at {csv_file}")
        return []
    
    return csv_file

def process_protein(
    config : TunnelingConfig,
    pdb_path: Path,
) -> list[Path]:
    
    prot_name = pdb_path.stem.removesuffix('_H')
    
    csv_file = get_predictions_csv(prot_name, config.pocket_predictions_dir)
    pocket_preds = pd.read_csv(csv_file, skipinitialspace=True)
    centers = pocket_preds[['center_x', 'center_y', 'center_z']]
    
    print(f"--- Processing {prot_name} ---")
    created_tunnels: list[Path] = []
    for pocket_idx in get_pocket_indices(config.pockets, pocket_preds):
        created_tunnels.extend(
            process_pocket(
                pocket_idx=pocket_idx,
                pocket_center_coords=centers.loc[pocket_idx],
                prot_name=prot_name,
                pdb_path=pdb_path,
                config=config,
            )
        )
    
    return created_tunnels

def run_tunneling(config: TunnelingConfig) -> list[Path]:
    pdb_dir = resolve_path(config.pdb_dir)
    output_dir = resolve_path(config.output_dir)
    caver_jar = resolve_path(config.caver_jar)


    if not pdb_dir.exists():
        raise FileNotFoundError(f"PDB directory not found: {pdb_dir}")
    if not caver_jar.exists():
        raise FileNotFoundError(f"CAVER jar not found: {caver_jar}")

    output_dir.mkdir(parents=True, exist_ok=True)
    save_tunnel_config(output_dir, config=config)

    created_tunnels: list[Path] = []
    for pdb_path in sorted(pdb_dir.glob("*.pdb")):
        created_tunnels.extend(
            process_protein(
                pdb_path=pdb_path,
                config=config
            )
        )
    return created_tunnels


def main() -> None:
    parser = create_parser()
    args = parser.parse_args()
    config = TunnelingConfig(
        pdb_dir=args.pdb_dir,
        pocket_predictions_dir=args.preds_dir,
        output_dir=args.output_dir,
        caver_jar=args.caver_jar,
        prep_receptor_dir=args.prep_receptor_dir,
        pockets=args.pockets,
        tunnels=args.tunnels,
        delta=args.delta,
        shell_depth=args.shell_depth,
        clustering_threshold=args.clustering_threshold,
        compute_tunnel_residues=args.compute_tunnel_residues,
        centroid_res_id=args.centroid_res_id,
        tol=args.tol
    )
    run_tunneling(config)


if __name__ == "__main__":
    main()
