import pandas as pd
import sys
import os
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method
import argparse

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--min-percent", type=float, required=True, help="Threshold in [0,1]")
    p.add_argument("--region", required=True, choices=["Sinushorn-SAnode-AtriumR"], help="Cardiac niche")
    # p.add_argument("--young-or-old", required=True, choices=["young","old"], help="Age group")
    return p.parse_args()

if __name__ == "__main__":
    args = parse_args()
    if not (0.0 <= args.min_percent <= 1.0):
        raise ValueError("--min-percent must be between 0 and 1")
    print(f"##### MIN_PERCENT={args.min_percent}, REGION={args.region} #####")
    # ... use args.min_percent / args.region / args.young_or_old ...

# Input files
cpdb_file_path = '/rfs/project/rfs-iCNyzSAaucw/kk837/databases/cellphonedb/v5.0.0/cellphonedb.zip'
meta_file_path = f'/rfs/project/rfs-iCNyzSAaucw/kk837/data_objects/Foetal/VisiumHD/Revision_Oct2025/cellphonedb/inputs/meta_{args.region}.txt'
counts_file_path = f'/rfs/project/rfs-iCNyzSAaucw/kk837/data_objects/Foetal/VisiumHD/Revision_Oct2025/cellphonedb/inputs/log_norm_counts_{args.region}.h5ad'
out_path = f'/rfs/project/rfs-iCNyzSAaucw/kk837/data_objects/Foetal/VisiumHD/Revision_Oct2025/cellphonedb/{args.region}'
if not os.path.exists(out_path):
    os.makedirs(out_path)

# run
cpdb_results = cpdb_statistical_analysis_method.call(
    cpdb_file_path = cpdb_file_path,                 # mandatory: CellphoneDB database zip file.
    meta_file_path = meta_file_path,                 # mandatory: tsv file defining barcodes to cell label.
    counts_file_path = counts_file_path,             # mandatory: normalized count matrix - a path to the counts file, or an in-memory AnnData object
    counts_data = 'hgnc_symbol',                     # defines the gene annotation in counts matrix.
    active_tfs_file_path = None,           # optional: defines cell types and their active TFs.
    microenvs_file_path = None,       # optional (default: None): defines cells per microenvironment.
    score_interactions = True,                       # optional: whether to score interactions or not. 
    iterations = 1000,                               # denotes the number of shufflings performed in the analysis.
    threshold = float(args.min_percent),                                 # defines the min % of cells expressing a gene for this to be employed in the analysis.
    threads = 5,                                     # number of threads to use in the analysis.
    debug_seed = 42,                                 # debug randome seed. To disable >=0.
    result_precision = 3,                            # Sets the rounding for the mean values in significan_means.
    pvalue = 0.1,                                   # P-value threshold to employ for significance.
    subsampling = False,                             # To enable subsampling the data (geometri sketching).
    subsampling_log = False,                         # (mandatory) enable subsampling log1p for non log-transformed data inputs.
    subsampling_num_pc = 100,                        # Number of componets to subsample via geometric skectching (dafault: 100).
    subsampling_num_cells = 1000,                    # Number of cells to subsample (integer) (default: 1/3 of the dataset).
    separator = '|',                                 # Sets the string to employ to separate cells in the results dataframes "cellA|CellB".
    debug = False,                                   # Saves all intermediate tables employed during the analysis in pkl format.
    output_path = out_path,                          # Path to save results.
    output_suffix = None                             # Replaces the timestamp in the output files by a user defined string in the  (default: None).
    )