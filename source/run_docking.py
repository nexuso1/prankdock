import argparse
import os
import sys
from glob import glob
from prody import *
from pathlib import Path
import subprocess

def dock_ligands(receptor_info : list[tuple[Path, Path]], ligands_folder, output_dir='../output', exhaustiveness=32):
    if not os.path.exists('../data/docking_files'):
        os.makedirs('../data/docking_files')

    if not os.path.exists('../output'):
        os.mkdir('../output')
    
    ligands = [Path(p).stem for p in glob(f'{ligands_folder}/*.pdbqt')]
    for receptor, config in receptor_info:
        receptor_id = receptor.stem.split('_')[0] # <id>_H_p<pocket_num>*
        pocket = receptor.stem.split('_')[2].removesuffix('.pdbqt')
        out_path = f'{output_dir}/{receptor_id}/{pocket}'.removesuffix('_prepared')
        if not os.path.exists(out_path):
            os.makedirs(out_path, exist_ok=True)
        elif all([os.path.exists(f'{out_path}/{lig}_out.pdbqt') for lig in ligands]):
            print(f'Already finished {receptor_id} pocket {pocket}, continuing..')
            continue

        print(f'Docking ligands to {receptor.stem}...')
        command = ' '.join([
            'vina',
            '--receptor', str(receptor),
            '--batch', f'{ligands_folder}/*.pdbqt',
            '--config', str(config),
            f'--exhaustiveness={exhaustiveness}',
            '--cpu', str(os.cpu_count()),
            '--dir', out_path
        ])
        subprocess.run(command, shell=True, check=True)
        print(f'Output saved to {out_path}')

def run_docking(args):
    receptors = list(glob(f'{args.dock_files_path}/**/*.pdbqt', recursive=True))
    receptors = [Path(r) for r in receptors]
    configs = [f'{rec.parent}/{rec.stem}.box.txt' for rec in receptors]
    receptor_info = list(zip(receptors, configs))
    dock_ligands(receptor_info, args.ligands_path, exhaustiveness=args.exhaustiveness)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-d', '--dock_files_path', default='../data/docking_files', help='Path to the folder where the preprared receptor .pdbqt files and vina configs are stored.')
    parser.add_argument('-l', '--ligands_path', default='../data/prepared_ligands', help='Folder containing prepared ligands for docking.')
    parser.add_argument('-e', '--exhaustiveness', type=int, default=32, help='Search exhaustiveness for Vina')
    parser.add_argument('-o', '--output_dir', default='../output')
    args = parser.parse_args()
    run_docking(args)
