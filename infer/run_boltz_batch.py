import yaml
import argparse
import subprocess
import re
import time
from pathlib import Path
from tqdm import tqdm

def parse_fasta(fasta_path):
    """Extract sequences from a FASTA file."""
    sequences = []
    current_id, current_seq = None, []
    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line: continue
            if line.startswith(">"):
                if current_id: sequences.append((current_id, "".join(current_seq)))
                current_id, current_seq = line[1:], []
            else:
                current_seq.append(line)
    if current_id: sequences.append((current_id, "".join(current_seq)))
    return sequences

def parse_smiles(smiles_path):
    """Parse 'ID: SMILES' from a file."""
    ligands = []
    if not smiles_path:
        return []
    with open(smiles_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or ":" not in line:
                continue
            lig_id, smiles = line.split(":", 1)
            ligands.append((lig_id.strip(), smiles.strip()))
    return ligands

def get_a3m_query(a3m_path):
    """Extract the first (query) sequence from an A3M file."""
    with open(a3m_path, 'r') as f:
        lines = f.readlines()
        if not lines: return None, None, [], 0
        header = lines[0].strip()
        seq_parts = []
        second_record_idx = len(lines)
        for i, line in enumerate(lines[1:], 1):
            if line.startswith(">"):
                second_record_idx = i
                break
            seq_parts.append(line.strip())
        return header, "".join(seq_parts), lines, second_record_idx

def ensure_msa_matches(fasta_seq, a3m_path, corrected_msa_dir):
    """Checks if MSA query matches FASTA; replaces it if not."""
    q_header, q_seq, all_lines, second_idx = get_a3m_query(a3m_path)
    
    # Remove any gaps or dots from MSA query for comparison
    clean_q_seq = q_seq.replace("-", "").replace(".", "")
    
    if clean_q_seq == fasta_seq:
        return a3m_path

    corrected_msa_dir.mkdir(parents=True, exist_ok=True)
    # Use a hash or the fasta_seq length in filename to avoid collisions
    new_a3m_path = corrected_msa_dir / f"corrected_{len(fasta_seq)}_{a3m_path.name}"
    
    if not new_a3m_path.exists():
        with open(new_a3m_path, 'w') as f:
            f.write(f"{all_lines[0].strip()}\n") # Original header
            f.write(f"{fasta_seq}\n")           # New sequence
            f.writelines(all_lines[second_idx:]) # Rest of the alignment
        
    return new_a3m_path

def main():
    parser = argparse.ArgumentParser(description="Automated Boltz Batch Prediction for FASTA variants and Ligands")
    parser.add_argument("--fasta_dir", required=True, help="Input FASTA directory")
    parser.add_argument("--a3m_dir", required=True, help="Input A3M directory")
    parser.add_argument("--smiles", help="Optional file containing 'ID: SMILES' for ligands")
    parser.add_argument("--predict_affinity", action="store_true", help="Enable affinity prediction for ligands")
    parser.add_argument("--affinity_samples", type=int, default=5, help="Number of diffusion samples for affinity prediction (default: 5)")
    parser.add_argument("--out_dir", required=True, help="Final output directory")
    parser.add_argument("--mamba_path", default="/home/ubuntu/miniforge3/bin/mamba", help="Path to mamba/conda")
    
    args = parser.parse_args()
    
    fasta_dir, a3m_dir = Path(args.fasta_dir), Path(args.a3m_dir)
    out_dir = Path(args.out_dir)
    config_dir = out_dir / "configs"
    predict_dir = out_dir / "predictions"
    corrected_msa_dir = out_dir / "corrected_msa"
    
    for d in [config_dir, predict_dir]: d.mkdir(parents=True, exist_ok=True)
    
    fasta_files = list(fasta_dir.glob("*.fasta")) + list(fasta_dir.glob("*.fa"))
    a3m_files = {f.stem: f for f in a3m_dir.glob("*.a3m")}
    ligands = parse_smiles(args.smiles)

    if not fasta_files:
        print(f"Error: No FASTA files found in {fasta_dir}")
        return

    universal_a3m = list(a3m_files.values())[0] if len(a3m_files) == 1 else None

    # Step 1: Pre-calculate total tasks (Cross-product if ligands are provided)
    all_tasks = []
    for f_path in fasta_files:
        a3m_match = a3m_files.get(f_path.stem) or universal_a3m
        if not a3m_match:
            continue
        variants = parse_fasta(f_path)
        for var_name, var_seq in variants:
            if ligands:
                for lig_id, smiles_str in ligands:
                    all_tasks.append({
                        "var_name": var_name,
                        "var_seq": var_seq,
                        "lig_id": lig_id,
                        "smiles": smiles_str,
                        "a3m_match": a3m_match,
                        "source_fasta": f_path.name
                    })
            else:
                all_tasks.append({
                    "var_name": var_name,
                    "var_seq": var_seq,
                    "lig_id": None,
                    "smiles": None,
                    "a3m_match": a3m_match,
                    "source_fasta": f_path.name
                })

    total_tasks = len(all_tasks)
    if total_tasks == 0:
        print("No valid tasks found.")
        return

    print(f"Found {total_tasks} prediction tasks across {len(fasta_files)} files.")

    # Step 2: Process with progress bar
    with tqdm(total=total_tasks, desc="Overall Progress", unit="task") as pbar:
        for task in all_tasks:
            var_name, var_seq = task["var_name"], task["var_seq"]
            lig_id, smiles = task["lig_id"], task["smiles"]
            a3m_match = task["a3m_match"]
            
            # Sanitization and task identification
            safe_var_name = re.sub(r'[^a-zA-Z0-9]', '_', var_name).strip('_')
            job_name = safe_var_name
            if lig_id:
                safe_lig_id = re.sub(r'[^a-zA-Z0-9]', '_', lig_id).strip('_')
                job_name = f"{safe_var_name}_{safe_lig_id}"
            
            pbar.set_description(f"Processing: {job_name}")
            
            # Ensure MSA matches for the protein chain
            msa_to_use = ensure_msa_matches(var_seq, a3m_match, corrected_msa_dir)

            # 1. Create YAML
            yaml_data = {
                "version": 1,
                "sequences": [
                    {"protein": {"id": "A", "sequence": var_seq, "msa": str(msa_to_use.resolve())}}
                ]
            }
            
            if smiles:
                yaml_data["sequences"].append({
                    "ligand": {"id": "B", "smiles": smiles}
                })
                
                if args.predict_affinity:
                    yaml_data["properties"] = [
                        {"affinity": {"binder": "B"}}
                    ]

            yaml_path = config_dir / f"{job_name}.yaml"
            with open(yaml_path, "w") as f:
                yaml.dump(yaml_data, f, sort_keys=False)

            # 2. Run Prediction
            cmd = [
                args.mamba_path, "run", "-n", "boltz", 
                "boltz", "predict", str(yaml_path),
                "--out_dir", str(predict_dir / job_name)
            ]
            
            if args.predict_affinity:
                cmd.extend(["--diffusion_samples_affinity", str(args.affinity_samples)])
            
            try:
                subprocess.run(cmd, check=True, capture_output=True, text=True)
            except subprocess.CalledProcessError as e:
                tqdm.write(f"!!! Error in {job_name}: {e.stderr}")
            
            pbar.update(1)

    print(f"\nBatch processing complete. Results are in: {out_dir.resolve()}")

if __name__ == "__main__":
    main()
