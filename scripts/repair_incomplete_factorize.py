#!/usr/bin/env python
"""
Repair script for cNMF incomplete factorization iterations.
This script detects incomplete factorization tasks and reruns them with skip_completed_runs=True.

Usage:
    python repair_incomplete_factorize.py --output-dir <dir> --name <name>
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
from pathlib import Path

def load_df_from_npz(filename):
    """Load dataframe from npz file"""
    try:
        with np.load(filename, allow_pickle=True) as f:
            obj = pd.DataFrame(**f)
        return obj
    except Exception as e:
        print(f"Error loading {filename}: {e}")
        return None

def update_nmf_iter_params(output_dir, name):
    """
    Update the nmf_replicate_parameters file to reflect actual completed status
    based on files that actually exist in the filesystem.
    """
    params_file = os.path.join(output_dir, name, 'cnmf_tmp', f'{name}.nmf_params.df.npz')
    
    if not os.path.exists(params_file):
        print(f"ERROR: Parameter file not found: {params_file}")
        return False
    
    # Load existing parameters
    replicate_params = load_df_from_npz(params_file)
    if replicate_params is None:
        return False
    
    # Update completed status based on actual files
    for i in replicate_params.index:
        k = replicate_params.at[i, 'n_components']
        iter_num = replicate_params.at[i, 'iter']
        spectra_file = os.path.join(output_dir, name, 'cnmf_tmp', 
                                   f'{name}.spectra.k_{k}.iter_{iter_num}.df.npz')
        
        if os.path.exists(spectra_file):
            replicate_params.at[i, 'completed'] = True
        else:
            replicate_params.at[i, 'completed'] = False
    
    # Save updated parameters
    try:
        np.savez_compressed(params_file, 
                          data=replicate_params.values, 
                          index=replicate_params.index.values, 
                          columns=replicate_params.columns.values)
        print(f"✓ Updated parameter file: {params_file}")
        return True
    except Exception as e:
        print(f"ERROR: Failed to save parameter file: {e}")
        return False

def detect_incomplete_tasks(output_dir, name):
    """
    Detect incomplete factorization tasks.
    Returns a summary of incomplete tasks grouped by k value.
    """
    params_file = os.path.join(output_dir, name, 'cnmf_tmp', f'{name}.nmf_params.df.npz')
    
    if not os.path.exists(params_file):
        print(f"ERROR: Parameter file not found: {params_file}")
        return None
    
    replicate_params = load_df_from_npz(params_file)
    if replicate_params is None:
        return None
    
    incomplete = replicate_params[replicate_params['completed'] == False]
    
    if len(incomplete) == 0:
        print("✓ All factorization tasks are complete!")
        return None
    
    print(f"\n{'='*70}")
    print(f"INCOMPLETE FACTORIZATION TASKS DETECTED: {len(incomplete)} tasks")
    print(f"{'='*70}")
    
    # Group by k value
    incomplete_by_k = {}
    for k in sorted(incomplete['n_components'].unique()):
        k_tasks = incomplete[incomplete['n_components'] == k]
        iters = sorted(k_tasks['iter'].values)
        incomplete_by_k[k] = iters
        print(f"\nk={k}: {len(k_tasks)} incomplete iterations")
        print(f"  Missing iterations: {iters}")
    
    print(f"\n{'='*70}\n")
    return incomplete_by_k

def repair_incomplete_factorize(output_dir, name, total_workers=1):
    """
    Repair incomplete factorization tasks by running factorize with skip_completed_runs=True.
    """
    print(f"\nStarting repair of incomplete factorization tasks...")
    print(f"Output dir: {output_dir}")
    print(f"Name: {name}")
    print(f"Total workers: {total_workers}")
    
    # Step 1: Update parameter file
    print(f"\nStep 1: Scanning filesystem to update task completion status...")
    if not update_nmf_iter_params(output_dir, name):
        print("ERROR: Failed to update parameter file")
        return False
    
    # Step 2: Detect incomplete tasks
    print(f"\nStep 2: Detecting incomplete tasks...")
    incomplete_by_k = detect_incomplete_tasks(output_dir, name)
    
    if incomplete_by_k is None:
        print("No incomplete tasks found. Exiting.")
        return True
    
    # Step 3: Run factorize with skip_completed_runs
    print(f"\nStep 3: Running factorize to repair incomplete iterations...")
    print(f"Using skip_completed_runs=True to skip already completed tasks")
    
    try:
        from cnmf import cNMF
        print(f"\nInitializing cNMF object...")
        cnmf_obj = cNMF(output_dir=output_dir, name=name)
        
        print(f"Running factorize with skip_completed_runs=True...")
        cnmf_obj.factorize(worker_i=0, total_workers=total_workers, skip_completed_runs=True)
        
        print(f"\n✓ Factorize repair completed!")
        
        # Verify completion
        print(f"\nStep 4: Verifying completion...")
        update_nmf_iter_params(output_dir, name)
        incomplete_after = detect_incomplete_tasks(output_dir, name)
        
        if incomplete_after is None:
            print("✓ All tasks are now complete!")
            return True
        else:
            print(f"⚠ Warning: Still {sum(len(v) for v in incomplete_after.values())} incomplete tasks")
            return False
            
    except ImportError:
        print(f"\nERROR: cNMF not installed. Please install it first.")
        return False
    except Exception as e:
        print(f"\nERROR during factorize: {e}")
        import traceback
        traceback.print_exc()
        return False

def main():
    parser = argparse.ArgumentParser(
        description='Repair incomplete cNMF factorization iterations',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Check for incomplete tasks and repair them
  python repair_incomplete_factorize.py --output-dir ./output --name my_analysis
  
  # Just detect without repairing (dry-run)
  python repair_incomplete_factorize.py --output-dir ./output --name my_analysis --detect-only
        """
    )
    
    parser.add_argument('--output-dir', required=True,
                       help='Output directory for cNMF analysis')
    parser.add_argument('--name', required=True,
                       help='Name of the cNMF analysis')
    parser.add_argument('--total-workers', type=int, default=1,
                       help='Total number of workers (default: 1)')
    parser.add_argument('--detect-only', action='store_true',
                       help='Only detect incomplete tasks, do not repair')
    
    args = parser.parse_args()
    
    # Step 1: Detect incomplete tasks
    if not update_nmf_iter_params(args.output_dir, args.name):
        sys.exit(1)
    
    incomplete_by_k = detect_incomplete_tasks(args.output_dir, args.name)
    
    if incomplete_by_k is None:
        print("\n✓ All tasks are complete. No repair needed.")
        sys.exit(0)
    
    if args.detect_only:
        print("\n(--detect-only flag set, not running repair)")
        sys.exit(0)
    
    # Step 2: Run repair
    success = repair_incomplete_factorize(args.output_dir, args.name, args.total_workers)
    
    sys.exit(0 if success else 1)

if __name__ == '__main__':
    main()

