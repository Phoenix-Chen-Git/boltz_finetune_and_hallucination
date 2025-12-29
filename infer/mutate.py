import argparse
from pathlib import Path

# Standard 20 amino acids
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"

def parse_fasta(fasta_path):
    """Simple FASTA parser."""
    with open(fasta_path, 'r') as f:
        lines = f.readlines()
        if not lines:
            raise ValueError("Empty FASTA file")
        header = lines[0].strip().lstrip(">")
        sequence = "".join(line.strip() for line in lines[1:])
    return header, sequence

def main():
    parser = argparse.ArgumentParser(description="Saturation Mutator for Protein Sequences")
    parser.add_argument("--fasta", required=True, help="Input FASTA file")
    parser.add_argument("--positions", required=True, help="Comma-separated 1-based indices (e.g. '10,25,100')")
    parser.add_argument("--out_dir", required=True, help="Output directory for the generated FASTA")
    
    args = parser.parse_args()
    
    # Parse inputs
    try:
        header, sequence = parse_fasta(args.fasta)
    except Exception as e:
        print(f"Error parsing FASTA: {e}")
        return
        
    seq_list = list(sequence)
    
    try:
        positions = [int(p.strip()) for p in args.positions.split(",")]
    except ValueError:
        print("Error: Positions must be a comma-separated list of integers.")
        return

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / f"{Path(args.fasta).stem}_saturated.fasta"

    variants_count = 0
    with open(out_file, 'w') as f:
        # Optionally include the wild-type first
        f.write(f">{header}_WT\n{sequence}\n")
        
        for pos in positions:
            if pos < 1 or pos > len(sequence):
                print(f"Warning: Position {pos} is out of range (1-{len(sequence)}). Skipping.")
                continue
            
            idx = pos - 1
            wt_aa = seq_list[idx]
            
            for mut_aa in AMINO_ACIDS:
                if mut_aa == wt_aa:
                    continue
                
                # Create mutation
                mut_seq = seq_list.copy()
                mut_seq[idx] = mut_aa
                
                # Format: OriginalHeader_Pos_WTtoMUT
                variant_name = f"{header}_{wt_aa}{pos}{mut_aa}"
                f.write(f">{variant_name}\n{''.join(mut_seq)}\n")
                variants_count += 1

    print(f"Success! Generated {variants_count} variants at {len(positions)} sites.")
    print(f"File saved to: {out_file}")

if __name__ == "__main__":
    main()

