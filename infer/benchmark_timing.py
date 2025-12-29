import subprocess
import time
import matplotlib.pyplot as plt
import os
import sys
from pathlib import Path

def run_boltz(yaml_path, out_dir, diff_samples=1, aff_samples=None):
    # Use 'boltz' directly as we'll run this script inside the mamba env
    cmd = [
        "boltz", "predict", str(yaml_path),
        "--out_dir", str(out_dir),
        "--diffusion_samples", str(diff_samples),
        "--override"
    ]
    if aff_samples is not None:
        cmd.extend(["--diffusion_samples_affinity", str(aff_samples)])
    
    print(f"Executing: {' '.join(cmd)}")
    start = time.time()
    result = subprocess.run(cmd, capture_output=True, text=True)
    end = time.time()
    
    if result.returncode != 0:
        print(f"Error running boltz:\nSTDOUT: {result.stdout}\nSTDERR: {result.stderr}")
        sys.exit(1)
        
    return end - start

def main():
    yaml_path = "/home/ubuntu/boltz_finetune_and_hallucination/infer/affinity_test/configs/variant_1_2_AG.yaml"
    base_out = Path("/home/ubuntu/boltz_finetune_and_hallucination/infer/benchmarks")
    base_out.mkdir(parents=True, exist_ok=True)

    print("Starting benchmarks...")
    
    # 1. Structure only, 1 sample
    print("\n[1/4] Running structure only (1 sample)...")
    t1 = run_boltz(yaml_path, base_out / "t1", diff_samples=1)
    
    # 2. Structure only, 5 samples
    print("\n[2/4] Running structure only (5 samples)...")
    t2 = run_boltz(yaml_path, base_out / "t2", diff_samples=5)
    
    # 3. Structure (1) + Affinity (5)
    print("\n[3/4] Running structure (1) + affinity (5 samples)...")
    t3 = run_boltz(yaml_path, base_out / "t3", diff_samples=1, aff_samples=5)
    
    # 4. Structure (1) + Affinity (10)
    print("\n[4/4] Running structure (1) + affinity (10 samples)...")
    t4 = run_boltz(yaml_path, base_out / "t4", diff_samples=1, aff_samples=10)

    # Derived timings
    struct_sample_time = (t2 - t1) / 4
    affinity_sample_time = (t4 - t3) / 5
    fixed_overhead = t1 - struct_sample_time

    print(f"\nResults:")
    print(f"Fixed Overhead (Preprocessing + Loading): {fixed_overhead:.2f}s")
    print(f"Time per Structure Sample: {struct_sample_time:.2f}s")
    print(f"Time per Affinity Sample: {affinity_sample_time:.2f}s")

    # Create a pie chart for a standard run: 1 structure sample, 5 affinity samples
    labels = ['Overhead', 'Structure Prediction (1 sample)', 'Affinity Prediction (5 samples)']
    sizes = [fixed_overhead, struct_sample_time, affinity_sample_time * 5]
    
    # Ensure no negative values due to noise
    sizes = [max(0.1, s) for s in sizes]

    plt.figure(figsize=(10, 7))
    plt.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=140)
    plt.title('Time Consumption Breakdown (1 Struct Sample, 5 Affinity Samples)')
    
    plot_path = "/home/ubuntu/boltz_finetune_and_hallucination/infer/timing_breakdown.png"
    plt.savefig(plot_path)
    print(f"\nPie plot saved to: {plot_path}")

if __name__ == "__main__":
    main()
