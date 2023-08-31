import os
import argparse
import logging
import sys


def adj_mtx_to_list(adj_mtx):
    """Converts adjacency matrix to list of edges.

    Parameters
    ----------
    adj_mtx : pd.DataFrame
        Adjacency matrix of the graph. Index and columns are the nodes.
    
    Returns
    -------
    pd.DataFrame
        List of edges with columns 'source', 'target', and 'weight'.
    """
    adj_mtx = adj_mtx.copy()
    adj_mtx = adj_mtx.stack().reset_index()
    adj_mtx.columns = ['source', 'target', 'weight']
    return adj_mtx


def main(args):

     # Logging setup
    log_format = "%(asctime)s - %(levelname)s - %(message)s"
    if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir)
    logging.basicConfig(filename=os.path.join(args.out_dir, "pyWGCNA_grn.log"), filemode='w', level=args.log.upper(), format=log_format)
    console = logging.StreamHandler()
    console.setLevel(args.log.upper())
    console.setFormatter(logging.Formatter(log_format))
    logging.getLogger().addHandler(console)

    # Parse rest of args
    h5ad_path = args.h5ad_path
    out_dir = args.out_dir
    tf_list = args.tf_list
    network_type = args.network_type
    layer = args.layer

    # Log args
    logging.info(f"Input h5ad file: {h5ad_path}")
    logging.info(f"Output directory: {out_dir}")
    logging.info(f"TF list: {tf_list}")
    logging.info(f"Network type: {network_type}")
    logging.info(f"Layer: {layer}")

    # Read in the h5ad and make sure the correct counts are used
    logging.info(f"Reading input file from {h5ad_path}")
    import scanpy as sc
    import pandas as pd
    adata = sc.read_h5ad(h5ad_path)
    if layer:
        logging.info(f"Using layer {layer} from adata")
        adata.X = adata.layers[layer]
        logging.info(f"Shape of adata: {adata.X.shape}")

    # Get the list of TFs
    logging.info(f"Reading TF list from {tf_list}")
    tfs = pd.read_csv(tf_list, header=None)[0].values

    # Get a dataframe for the expression data
    dat = adata.to_df()
    logging.info("Converted adata to dataframe.")
    logging.info(f"Shape of dataframe: {dat.shape}")
    
    # Run the PyWGCNA analysis
    logging.info("Starting a PyWGCNA GRN run")
    import PyWGCNA
    
    # Pick a soft threshold automatically
    logging.info("Picking soft threshold...")
    power, sft = PyWGCNA.WGCNA.pickSoftThreshold(data=dat, networkType=network_type)
    sft.to_csv(os.path.join(out_dir, "pick_soft_threshold.tsv"), sep="\t")
    logging.info(f"Soft threshold: {power}")

    # Calculate an adjacency matrix
    logging.info("Calculating adjacency matrix...")
    adjacency = PyWGCNA.WGCNA.adjacency(dat, power=power, adjacencyType=network_type)

    # Get an adjacency list for this and save it
    logging.info("Saving adjacency list...")
    adjacency_df = pd.DataFrame(adjacency, columns=dat.columns.values, index=dat.columns.values)
    corr_adj_list = adj_mtx_to_list(adjacency_df)
    corr_adj_list = corr_adj_list[corr_adj_list['weight'] != 0]
    corr_adj_list = corr_adj_list[corr_adj_list["source"].isin(tfs)]
    corr_adj_list['minGene'] = corr_adj_list[['source', 'target']].min(axis=1)
    corr_adj_list['maxGene'] = corr_adj_list[['source', 'target']].max(axis=1)
    corr_adj_list = corr_adj_list.drop_duplicates(subset=['minGene', 'maxGene'])
    corr_adj_list = corr_adj_list.drop(columns=['minGene', 'maxGene'])
    corr_adj_list.to_csv(os.path.join(out_dir, "corr_adj_list.tsv"), index=False, sep="\t")
    logging.info(f"Saved {corr_adj_list.shape[0]} edges to {os.path.join(out_dir, 'corr_adj_list.tsv')}")

    # Get the topological overlap matrix
    logging.info("Calculating TOM...")
    TOM = PyWGCNA.WGCNA.TOMsimilarity(adjacency, TOMType=network_type)
    
    # Clean up the TOM convert to an adjacency list and save it
    logging.info("Saving TOM adjacency list...")
    TOM.columns = dat.columns.values
    TOM.index = dat.columns.values
    tom_adj_list = adj_mtx_to_list(TOM)
    tom_adj_list = tom_adj_list[tom_adj_list['weight'] != 0]
    tom_adj_list = tom_adj_list[tom_adj_list["source"].isin(tfs)]
    tom_adj_list['minGene'] = tom_adj_list[['source', 'target']].min(axis=1)
    tom_adj_list['maxGene'] = tom_adj_list[['source', 'target']].max(axis=1)
    tom_adj_list = tom_adj_list.drop_duplicates(subset=['minGene', 'maxGene'])
    tom_adj_list = tom_adj_list.drop(columns=['minGene', 'maxGene'])
    tom_adj_list.to_csv(os.path.join(out_dir, "tom_adj_list.tsv"), index=False, sep="\t")
    logging.info(f"Saved {tom_adj_list.shape[0]} edges to {os.path.join(out_dir, 'tom_adj_list.tsv')}")
    
    # Done!
    logging.info("Done!")


if __name__ == "__main__":
    """Parse arguments and handle the main execution."""
    parser = argparse.ArgumentParser(description="Process the dataset and generate adjacency lists.")
    parser.add_argument("--h5ad_path", type=str, required=True, help="Path to the input h5ad file")
    parser.add_argument("--out_dir", type=str, required=True, help="Directory to save the results")
    parser.add_argument("--tf_list", type=str, required=True, help="Path to the list of TFs")
    parser.add_argument("--network_type", type=str, default="signed", choices=["signed", "unsigned"], help="Type of network for the analysis")
    parser.add_argument("--layer", type=str, default=None, help="Layer of the adata to use, if any")
    parser.add_argument("--log", type=str, default="info", choices=["debug", "info", "warning", "error", "critical"], help="Set the logging level")

    # Run with args
    args = parser.parse_args()
    main(args)
