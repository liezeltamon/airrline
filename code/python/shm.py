import pandas as pd
import numpy as np

def calculate_mutation_ratios(df):
    """
    Calculate clinically relevant CDR/FWR mutation ratios from immunoglobulin repertoire data.
    
    Parameters:
    -----------
    df : pd.DataFrame
        DataFrame containing mutation count/frequency columns with pattern:
        mu_{metric}_{region}_{mutation_type}_{isotype}
        
    Returns:
    --------
    pd.DataFrame
        Original dataframe with added ratio columns
    """
    
    def safe_ratio(numerator, denominator):
        """Calculate ratio with safe handling of division by zero"""
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = numerator / denominator
            ratio = np.where(np.isinf(ratio) | np.isnan(ratio), np.nan, ratio)
        return ratio
    
    # Create copy to avoid modifying original
    result_df = df.copy()
    
    # Check which metric types are available
    mu_cols = [col for col in df.columns if col.startswith("mu_count_") or col.startswith("mu_freq_")]
    has_count = any(col.startswith("mu_count_") for col in df.columns)
    has_freq = any(col.startswith("mu_freq_") for col in df.columns)
    
    def calc_ratios_for_type(metric_type):
        """Calculate ratios for either 'count' or 'freq' metrics"""
        
        # Overall CDR vs FWR ratios
        # Total CDR (CDR1 + CDR2)
        cdr_cols = [col for col in df.columns if 
                    col.startswith(f"mu_{metric_type}_cdr") and 
                    col.endswith(("_r", "_s")) and
                    not any(iso in col for iso in ["_IGH", "_IGK", "_IGL"])]
        cdr_total = df[cdr_cols].sum(axis=1) if cdr_cols else pd.Series(0, index=df.index)
        
        # Total FWR (FWR1 + FWR2 + FWR3)
        fwr_cols = [col for col in df.columns if 
                    col.startswith(f"mu_{metric_type}_fwr") and 
                    col.endswith(("_r", "_s")) and
                    not any(iso in col for iso in ["_IGH", "_IGK", "_IGL"])]
        fwr_total = df[fwr_cols].sum(axis=1) if fwr_cols else pd.Series(0, index=df.index)
        
        result_df[f"ratio_{metric_type}_cdr_fwr"] = safe_ratio(cdr_total, fwr_total)
        
        # R vs S mutation ratios
        # CDR R vs S
        cdr_r_cols = [col for col in df.columns if 
                      col.startswith(f"mu_{metric_type}_cdr") and 
                      col.endswith("_r") and
                      not any(iso in col for iso in ["_IGH", "_IGK", "_IGL"])]
        cdr_r = df[cdr_r_cols].sum(axis=1) if cdr_r_cols else pd.Series(0, index=df.index)
        
        cdr_s_cols = [col for col in df.columns if 
                      col.startswith(f"mu_{metric_type}_cdr") and 
                      col.endswith("_s") and
                      not any(iso in col for iso in ["_IGH", "_IGK", "_IGL"])]
        cdr_s = df[cdr_s_cols].sum(axis=1) if cdr_s_cols else pd.Series(0, index=df.index)
        
        result_df[f"ratio_{metric_type}_cdr_r_s"] = safe_ratio(cdr_r, cdr_s)
        
        # FWR R vs S
        fwr_r_cols = [col for col in df.columns if 
                      col.startswith(f"mu_{metric_type}_fwr") and 
                      col.endswith("_r") and
                      not any(iso in col for iso in ["_IGH", "_IGK", "_IGL"])]
        fwr_r = df[fwr_r_cols].sum(axis=1) if fwr_r_cols else pd.Series(0, index=df.index)
        
        fwr_s_cols = [col for col in df.columns if 
                      col.startswith(f"mu_{metric_type}_fwr") and 
                      col.endswith("_s") and
                      not any(iso in col for iso in ["_IGH", "_IGK", "_IGL"])]
        fwr_s = df[fwr_s_cols].sum(axis=1) if fwr_s_cols else pd.Series(0, index=df.index)
        
        result_df[f"ratio_{metric_type}_fwr_r_s"] = safe_ratio(fwr_r, fwr_s)
        
        # By isotype (IGH, IGK, IGL)
        for isotype in ["IGH", "IGK", "IGL"]:
            # CDR vs FWR by isotype
            cdr_iso_cols = [col for col in df.columns if 
                           col.startswith(f"mu_{metric_type}_cdr") and 
                           col.endswith(f"_{isotype}")]
            cdr_iso = df[cdr_iso_cols].sum(axis=1) if cdr_iso_cols else pd.Series(0, index=df.index)
            
            fwr_iso_cols = [col for col in df.columns if 
                           col.startswith(f"mu_{metric_type}_fwr") and 
                           col.endswith(f"_{isotype}")]
            fwr_iso = df[fwr_iso_cols].sum(axis=1) if fwr_iso_cols else pd.Series(0, index=df.index)
            
            result_df[f"ratio_{metric_type}_cdr_fwr_{isotype}"] = safe_ratio(cdr_iso, fwr_iso)
            
            # R vs S in CDR by isotype
            cdr_r_iso_cols = [col for col in df.columns if 
                             col.startswith(f"mu_{metric_type}_cdr") and 
                             col.endswith(f"_r_{isotype}")]
            cdr_r_iso = df[cdr_r_iso_cols].sum(axis=1) if cdr_r_iso_cols else pd.Series(0, index=df.index)
            
            cdr_s_iso_cols = [col for col in df.columns if 
                             col.startswith(f"mu_{metric_type}_cdr") and 
                             col.endswith(f"_s_{isotype}")]
            cdr_s_iso = df[cdr_s_iso_cols].sum(axis=1) if cdr_s_iso_cols else pd.Series(0, index=df.index)
            
            result_df[f"ratio_{metric_type}_cdr_r_s_{isotype}"] = safe_ratio(cdr_r_iso, cdr_s_iso)
            
            # R vs S in FWR by isotype
            fwr_r_iso_cols = [col for col in df.columns if 
                             col.startswith(f"mu_{metric_type}_fwr") and 
                             col.endswith(f"_r_{isotype}")]
            fwr_r_iso = df[fwr_r_iso_cols].sum(axis=1) if fwr_r_iso_cols else pd.Series(0, index=df.index)
            
            fwr_s_iso_cols = [col for col in df.columns if 
                             col.startswith(f"mu_{metric_type}_fwr") and 
                             col.endswith(f"_s_{isotype}")]
            fwr_s_iso = df[fwr_s_iso_cols].sum(axis=1) if fwr_s_iso_cols else pd.Series(0, index=df.index)
            
            result_df[f"ratio_{metric_type}_fwr_r_s_{isotype}"] = safe_ratio(fwr_r_iso, fwr_s_iso)
        
        # Individual region ratios
        # CDR1 vs FWR1, CDR2 vs FWR2
        for region_num in [1, 2]:
            cdr_region_cols = [col for col in df.columns if 
                              col.startswith(f"mu_{metric_type}_cdr{region_num}_") and
                              not any(iso in col for iso in ["_IGH", "_IGK", "_IGL"])]
            cdr_region = df[cdr_region_cols].sum(axis=1) if cdr_region_cols else pd.Series(0, index=df.index)
            
            fwr_region_cols = [col for col in df.columns if 
                              col.startswith(f"mu_{metric_type}_fwr{region_num}_") and
                              not any(iso in col for iso in ["_IGH", "_IGK", "_IGL"])]
            fwr_region = df[fwr_region_cols].sum(axis=1) if fwr_region_cols else pd.Series(0, index=df.index)
            
            result_df[f"ratio_{metric_type}_cdr{region_num}_fwr{region_num}"] = safe_ratio(cdr_region, fwr_region)
    
    # Calculate for counts if available
    if has_count:
        calc_ratios_for_type("count")
    
    # Calculate for frequencies if available
    if has_freq:
        calc_ratios_for_type("freq")
    
    return result_df

# Usage example:
# df_with_ratios = calculate_mutation_ratios(df)
# print(df_with_ratios.filter(regex='^ratio_').columns)

# %%
# import os
# import pandas as pd
# os.chdir("/gpfs3/well/immune-rep/users/yfg436/git/sle")
# mutations_path = "results/airr/quantify_mutations/ig_rearranged/all/mutations.csv"
# df = pd.read_csv(mutations_path, index_col=0)
# df_with_ratios = calculate_mutation_ratios(df)
# %%
