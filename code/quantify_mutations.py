# %% Quantify mutations in BCR chains
# env: dandelion
# sbatch -J quantify_mutations -p long,short --output=%x.log.out --error=%x.log.err --wrap='python quantify_mutations.py'

import os
os.chdir("../..")
import dandelion as ddl
import pandas as pd
from functools import reduce, partial

# %% Parameters

contig_dir = "code/airr/reannotate/contigs_separate_directories"
ignore_files = ["meta.csv", "reannotate_logs", "SLE_batch_2_pool_1", "SLE_batch_2_pool_2"]
library_type_id = "ig_rearranged"
ddl_prefix = "all"

out_dir = os.path.join(
    "results", "airr", "quantify_mutations", library_type_id, ddl_prefix
)
os.makedirs(out_dir, exist_ok=True)

# %% Function

# Merge data frame in a list on a common column
def merge_dfs(left, right):
    return pd.merge(
        left, right, left_index=True, right_index=True, suffixes=("", "_y"), how="outer"
    )
merge_dfs_partial = partial(merge_dfs)

# %% ----- MAIN -----

# %% Concatenate all ig contigs in one dandelion object

src_ids = list(set([id for id in os.listdir(os.path.join(contig_dir, library_type_id, ddl_prefix))]) - set(ignore_files))
assert len(src_ids) == len(set(src_ids)), "Duplicate src ids found"

ddl_list = []
for src_id in src_ids:
    print(f"Loading {src_id} ...")
    vdj_ddl = ddl.Dandelion(os.path.join(contig_dir, library_type_id, ddl_prefix, src_id, src_id, "dandelion", "all_contig_dandelion.tsv"))
    ddl_list.append(vdj_ddl)
vdj = ddl.concat(ddl_list)

assert vdj.metadata.index.is_unique, "Duplicate cell barcodes found in the concatenated dandelion object"

# %% Quantify mutations

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

# %% Save

mu_df.reset_index(names="index_unique", inplace=True)
mu_df.to_csv(os.path.join(out_dir, "mutations.csv"), index=False)

# mu_df.to_hdf(
#     os.path.join(out_dir, "mutations.h5"),
#     key="x"
# )
