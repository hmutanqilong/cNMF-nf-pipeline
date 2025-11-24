#!/usr/bin/env python
"""
Check if all expected factorize iterations are complete.
Exit codes:
  0 = All complete
  1 = Incomplete (can be repaired)
  2 = Error
"""

import os
import sys
import argparse

def check_completeness(output_dir, name, expected_iters, k_values):
    """
    Check if all (k, iter) combinations have corresponding files.
    Returns:
      0 if all complete
      1 if incomplete
      2 if error
    """
    missing_count = 0
    missing_by_k = {}
    
    for k in k_values:
        missing_iters = []
        for iter_num in range(expected_iters):
            spectra_file = os.path.join(
                output_dir, name, 'cnmf_tmp',
                f'{name}.spectra.k_{k}.iter_{iter_num}.df.npz'
            )
            if not os.path.exists(spectra_file):
                missing_iters.append(iter_num)
                missing_count += 1
        
        if missing_iters:
            missing_by_k[k] = missing_iters
    
    if missing_count == 0:
        print(f"✓ All {len(k_values)} k-values have all {expected_iters} iterations")
        return 0
    
    print(f"\n{'='*70}")
    print(f"❌ INCOMPLETE FACTORIZATIONS DETECTED: {missing_count} missing files")
    print(f"{'='*70}")
    for k in sorted(missing_by_k.keys()):
        print(f"\nk={k}: Missing {len(missing_by_k[k])} iterations")
        print(f"  {missing_by_k[k]}")
    print(f"\n{'='*70}\n")
    
    return 1

def main():
    parser = argparse.ArgumentParser(
        description='Check factorization completeness'
    )
    parser.add_argument('--output-dir', required=True)
    parser.add_argument('--name', required=True)
    parser.add_argument('--expected-iters', type=int, required=True)
    parser.add_argument('--k-values', required=True)
    
    args = parser.parse_args()
    k_values = [int(k) for k in args.k_values.split(',')]
    
    result = check_completeness(
        args.output_dir,
        args.name,
        args.expected_iters,
        k_values
    )
    sys.exit(result)

if __name__ == '__main__':
    main()
