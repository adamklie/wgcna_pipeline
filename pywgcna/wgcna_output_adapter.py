import os
import glob
import logging
import argparse
import time

logging.basicConfig(level=logging.INFO)

def main(args):
    start_time = time.time()
    logging.info("Adapter started.")
    import pandas as pd
    import numpy as np

    # Read in
    corr_file = glob.glob(os.path.join(args.wgcna_out_dir, "corr*adj_list.tsv"))[0]
    tom_file = glob.glob(os.path.join(args.wgcna_out_dir, "tom*adj_list.tsv"))[0]
    corr_adj = pd.read_csv(corr_file, sep="\t", low_memory=False)
    tom_adj = pd.read_csv(tom_file, sep="\t", low_memory=False)
    
    # Clean up
    corr_adj.rename(columns={"weight": "weight_signed"}, inplace=True)
    tom_adj.rename(columns={"weight": "weight_unsigned"}, inplace=True)
    corr_adj["weight_unsigned"] = corr_adj["weight_signed"].abs()
    tom_adj["weight_signed"] = np.nan
    corr_adj["weight_minmax_normalized"] = (corr_adj["weight_unsigned"] - corr_adj["weight_unsigned"].min()) / (corr_adj["weight_unsigned"].max() - corr_adj["weight_unsigned"].min())
    tom_adj["weight_minmax_normalized"] = (tom_adj["weight_unsigned"] - tom_adj["weight_unsigned"].min()) / (tom_adj["weight_unsigned"].max() - tom_adj["weight_unsigned"].min())
    corr_adj["p"] = np.nan
    corr_adj["-logp"] = np.nan
    corr_adj["description"] = np.nan
    tom_adj["p"] = np.nan
    tom_adj["-logp"] = np.nan
    tom_adj["description"] = np.nan

    # Save
    corr_adj.to_csv(os.path.join(args.wgcna_out_dir, f"corr{args.out_name}"), sep="\t", index=False)
    tom_adj.to_csv(os.path.join(args.wgcna_out_dir, f"tom{args.out_name}"), sep="\t", index=False)
    end_time = time.time()
    logging.info(f"Process completed. Results saved at {args.wgcna_out_dir}. Duration: {end_time - start_time:.2f} seconds")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Adapt WGCNA adjacencies.")
    parser.add_argument('--wgcna_out_dir', type=str, required=True, help='Path to the output directory.')
    parser.add_argument('--out_name', type=str, default="_net.tsv", help='Name of the output file.')
    
    args = parser.parse_args()
    main(args)
