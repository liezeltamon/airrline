# %% Quantify mutations in BCR chains
# env: dandelion
# srun -p long,short --output=log.out --error=log.err python quantify_mutations.py --library_type=ig &

import os
import argparse
import dandelion as ddl
import pandas as pd
from functools import reduce, partial

os.chdir("../..")

# %% Parameters

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Quantify mutations.")
parser.add_argument("--library_type", type=str, required=True, help="ddl.pp.check_contigs(..., library_type).")
args = parser.parse_args()

library_type = args.library_type
#library_type = "ig"

airr_path = os.path.join("results/airr/postprocess/filtered", library_type, "airr_rearrangement.tsv")
out_dir = os.path.join("results", "airr", "quantify_mutations", library_type)
os.makedirs(out_dir, exist_ok=True)

# %% Function

# Merge data frame in a list on a common column
def merge_dfs(left, right):
    return pd.merge(
        left, right, left_index=True, right_index=True, suffixes=("", "_y"), how="outer"
    )
merge_dfs_partial = partial(merge_dfs)

# %% ----- MAIN -----

# Load reannotated AIRR to adata
vdj = ddl.read_airr(airr_path)

mu_df_list = []
for bool in [True, False]:

    ddl.pp.quantify_mutations(
        vdj,
        # "Whether to return the results for heavy chain and light chain separately"
        split_locus=True,
        frequency = bool,
        # "Whether to return the results for replacement and silent mutations separately"
        combine=False,
        # For region definitions, see https://shazam.readthedocs.io/en/stable/topics/IMGT_SCHEMES/
        region_definition="IMGT_V_BY_REGIONS"
    )
    mu_df_list.append(vdj.metadata.filter(like="mu_").copy())

for bool in [True, False]:

    ddl.pp.quantify_mutations(
        vdj,
        # "Whether to return the results for heavy chain and light chain separately"
        split_locus=False,
        frequency = bool,
        # "Whether to return the results for replacement and silent mutations separately"
        combine=True
    )
    mu_df_list.append(vdj.metadata.filter(like="mu_").copy())

mu_df = reduce(merge_dfs_partial, mu_df_list)
# Drop the duplicated columns with the "_y" suffix
mu_df = mu_df.loc[:, ~mu_df.columns.str.endswith("_y")]
mu_df.sort_index(axis=1, inplace=True)

mu_df.to_hdf(
    os.path.join(out_dir, "mutations.h5"),
    key="x"
)
