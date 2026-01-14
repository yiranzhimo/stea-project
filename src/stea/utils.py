import os
import scanpy as sc
import pandas as pd
import numpy as np
import seaborn as sns
from scipy.stats import rankdata
from scipy.sparse import issparse
from matplotlib import pyplot as plt 
from sklearn.metrics.pairwise import cosine_similarity 
from scipy.special import kl_div
from statsmodels.stats.multitest import multipletests


def qcplot(adata, folder="01.QC"):
    """
    Parameters:
    -----------
    adata : AnnData
        Spatial transcriptomics data object containing gene expression matrix
        - adata.X: Gene expression matrix (cells x genes)
        - adata.var_names: Gene name list
        - adata.obs_names: Cell name list
    
    folder : str, default="01.QC"
        Output folder name for saving QC plots and results
    
    Returns:
    --------
    None
        The function does not return a value, but it will:
        1. Detect and report data normalization status
        2. Add QC metrics columns to adata.obs
        3. Generate QC visualization plots and save to the specified folder
    """
    if not os.path.exists(folder):
        os.makedirs(folder)
        print(f"Folder '{folder}' was created.")
    else:
        print(f"Folder '{folder}' already exists.")

    adata.var_names_make_unique()
    adata.var["mt"] = adata.var_names.str.startswith("^MT-|^mt-|^Mt-")
    # Detect data normalization status
    if issparse(adata.X):
        X = adata.X.toarray()
    else:
        X = adata.X
    has_decimals = np.any(X % 1 != 0)  # Check if there are decimals
    has_negatives = np.any(X < 0)      # Check if there are negative values
    if has_negatives:
        print("Data has been scaled! Please use raw data or nomalized data for analysis")
    elif has_decimals:
        print("Data has been normalized! ")
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
        fig, axs = plt.subplots(1, 2, figsize=(8, 4))
        sns.histplot(adata.obs["total_counts"], kde=False,bins=50, ax=axs[0])
        sns.histplot(adata.obs["n_genes_by_counts"], kde=False, bins=50, ax=axs[1])
        fig.savefig(folder+"/QC.png")
    else:
        print("Data has not been normalized! ")
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
    fig, axs = plt.subplots(1, 2, figsize=(8, 4))
    sns.histplot(adata.obs["total_counts"], kde=False,bins=50, ax=axs[0])
    sns.histplot(adata.obs["n_genes_by_counts"], kde=False, bins=50, ax=axs[1])
    fig.savefig(folder+"/QC.png")


def marker_rank_analysis(adata, genesets, folder="02.MarkerRank"):
    """
    Calculate marker gene expression ranks for each cell type in gene sets and generate visualization plots.
    
    Parameters:
    -----------
    adata : AnnData
        Single-cell data object containing gene expression matrix
        - adata.X: Gene expression matrix (cells x genes)
        - adata.var_names: Gene names
    
    genesets : dict
        Gene set dictionary with format {cell_type: [marker_gene_list]}
        Example: {
            'T_cell': ['CD3D', 'CD3E', 'CD8A', ...],
            'B_cell': ['CD19', 'MS4A1', 'CD79A', ...],
            'NK_cell': ['GNLY', 'NKG7', 'GZMB', ...]
        }
    
    folder : str, default="02.MarkerRank"
        Output folder name for saving ranking results and plots
    
    Returns:
    --------
    dict
        Dictionary containing marker gene ranking results for each cell type
        Format: {cell_type: pd.DataFrame}
    """
    if not os.path.exists(folder):
        os.makedirs(folder)
        print(f"Folder '{folder}' was created.")
    else:
        print(f"Folder '{folder}' already exists.")
    
    # Get cell type list (from genesets keys)
    celltypes = list(genesets.keys())
    print(f"Found {len(celltypes)} cell types from genesets: {celltypes}")
    
    # Save original count data
    if 'counts' not in adata.layers:
        adata.layers['counts'] = adata.X.copy()
        print("Saved original counts to adata.layers['counts']")
    
    # Convert to DataFrame and transpose (genes x cells)
    data = adata.to_df().T
    print(f"Data shape: {data.shape} (genes x cells)")
    
    # Rank gene expression for each cell (descending, lower rank = higher expression)
    data_ranked = data.rank(ascending=False)
    print("Calculated expression ranks for each cell")
    
    # Store all results
    all_results = {}
    
    # Analyze marker genes for each cell type
    for celltype in genesets.keys():
        print(f"\nAnalyzing marker genes for {celltype}...")
        
        # Get marker genes for this cell type
        marker_genes = genesets[celltype]
        
        # Check if marker genes exist in the data
        available_genes = [gene for gene in marker_genes if gene in data_ranked.index]
        missing_genes = [gene for gene in marker_genes if gene not in data_ranked.index]
        
        if missing_genes:
            print(f"Warning: {len(missing_genes)} genes not found in data: {missing_genes}")
        
        if not available_genes:
            print(f"Error: No marker genes found for {celltype}")
            continue
            
        print(f"Using {len(available_genes)} available marker genes")
        
        # Extract ranking data for marker genes
        sub_data = data_ranked.loc[available_genes]
        
        # Calculate median rank for each gene
        sub_data['median_rank'] = sub_data.median(axis=1)
        
        # Sort by median rank (ascending, lower rank is better)
        sorted_data = sub_data.sort_values('median_rank', ascending=True)
        
        # Save ranking results
        rank_results = sorted_data.drop('median_rank', axis=1).copy()
        rank_results['median_rank'] = sorted_data['median_rank']
        rank_results['gene'] = rank_results.index
        
        # Save to CSV
        rank_results.to_csv(f"{folder}/{celltype}_marker_ranks.csv")
        all_results[celltype] = rank_results
        
        # Reshape data for plotting
        plot_data = sorted_data.drop('median_rank', axis=1).copy()
        plot_data["gene"] = plot_data.index
        plot_data_melted = pd.melt(plot_data, id_vars="gene", var_name="cell", value_name="rank")
        
        # Create visualization plot
        plt.figure(figsize=(12, max(8, len(available_genes) * 0.3)))
        sns.stripplot(x="rank", y="gene", data=plot_data_melted, orient='h', size=1, alpha=0.6)
        
        # Add reference lines
        reference_lines = [100, 500, 1000, 5000, 10000]
        for line in reference_lines:
            if line < data_ranked.shape[1]:  # Ensure reference line is within data range
                plt.axvline(x=line, color='black', linestyle='--', alpha=0.7)
        
        plt.xlabel('Expression Rank (lower = higher expression)')
        plt.ylabel('Marker Genes')
        plt.title(f'{celltype} Marker Gene Expression Ranks')
        plt.tight_layout()
        
        # Save plot
        plt.savefig(f"{folder}/{celltype}_marker_ranks.png", dpi=300, bbox_inches='tight')
        plt.close()
        
    
    # Create summary statistics
    summary_stats = []
    for celltype, results in all_results.items():
        median_ranks = results['median_rank']
        summary_stats.append({
            'celltype': celltype,
            'n_markers': len(results),
            'min_rank': median_ranks.min(),
            'max_rank': median_ranks.max(),
            'mean_rank': median_ranks.mean(),
            'median_rank': median_ranks.median(),
            'top_gene': results.loc[results['median_rank'].idxmin(), 'gene']
        })
    return all_results

def stgsea(adata, genesets, alpha=0.25, norm=False, folder="03.ES"):
    """
    Parameters:
    -----------
    adata : AnnData
        Single-cell data object containing gene expression matrix
        - adata.X: Gene expression matrix (cells x genes)
        - adata.var_names: Gene names
        - adata.obs_names: Cell names
    
    genesets : dict
        Gene set dictionary with format {gene_set_name: [gene_list]}
        Example: {
            'Cell_Cycle': ['CCNA2', 'CCNB1', 'CDK1', ...],
            'Immune_Response': ['IFNG', 'TNF', 'IL1B', ...],
            'Apoptosis': ['BAX', 'BCL2', 'CASP3', ...]
        }
    
    alpha : float, default=0.25
        Weight parameter in enrichment score calculation
        - Smaller alpha values (e.g., 0.25) focus more on highly expressed genes
        - Larger alpha values (e.g., 1.0) give equal weight to all genes
        - Recommended range: 0.1-1.0
    
    norm : bool, default=False
        Whether to normalize enrichment scores
        - True: Normalize enrichment scores to 0-1 range
        - False: Keep original enrichment scores
    
    folder : str, default="03.ES"
        Output folder name for saving enrichment score results
    
    Returns:
    --------
    pd.DataFrame
        Enrichment score matrix with shape (number of gene sets x number of cells)
        - Row index: Gene set names
        - Column index: Cell names
        - Values: Enrichment scores
    
    Notes:
    ------
    1. Enrichment score calculation is based on gene expression ranking statistics
    2. Positive values indicate gene set enrichment (high expression) in that cell
    3. Negative values indicate gene set depletion (low expression) in that cell
    4. Larger absolute values indicate stronger enrichment
    """
    if not os.path.exists(folder):
        os.makedirs(folder)
        print(f"Folder '{folder}' was created.")
    else:
        print(f"Folder '{folder}' already exists.")
    row_names = adata.var_names
    num_genes = adata.shape[1]
    gene_sets = [np.where(row_names.isin(genesets[key]))[0] for key in genesets.keys()]
    if issparse(adata.X):
        dense_X = adata.X.toarray()
    else:
        dense_X= adata.X
    R = rankdata(dense_X, axis=1)
    R = 10000 * R/R.shape[1]
    es = np.array([
        np.apply_along_axis(
            lambda R_col: np.array([
                np.sum((step_cdf_pos - step_cdf_neg) / num_genes)
                for gene_ranks in [np.argsort(-R_col)]
                for indicator_pos in [np.isin(gene_ranks, gene_set)]
                for indicator_neg in [~indicator_pos]
                for step_cdf_pos, step_cdf_neg in [
                    (
                        np.cumsum((R_col[gene_ranks] * indicator_pos) ** alpha) / np.sum((R_col[gene_ranks] * indicator_pos) ** alpha),
                        np.cumsum(indicator_neg) / np.sum(indicator_neg)
                    )
                ]
            ]), axis=1, arr=R)
        for gene_set in gene_sets
    ])

    es = np.squeeze(es)

    if len(gene_sets) == 1:
        es = es.reshape(1, -1)

    if norm:
        es = es / (np.nanmax(es) - np.nanmin(es))
    es = pd.DataFrame(es, index=list(genesets.keys()), columns=adata.obs_names)
    
    es.to_csv(folder+"/stgsea.csv")
    max_rows = es.idxmax().to_frame()
    max_rows.columns = ['celltype']
    max_rows.to_csv(folder+"/celltype_prediction.csv")

    return max_rows


def calculate_gene_cosine_similarity(adata, target_gene, top_n=None, folder="04.CosineSimilarity"):
    """
    Calculate cosine similarity between a specified gene and all other genes.
    
    Parameters:
    -----------
    adata : AnnData
        Spatial transcriptomics data object containing gene expression matrix
        - adata.X: Gene expression matrix (cells x genes)
        - adata.var_names: Gene name list
        - adata.obs_names: Cell name list
    target_gene : str
        Target gene name
    top_n : int, optional, default=None
        Return top n genes with highest similarity. If None, return all genes
    folder : str, default="04.CosineSimilarity"
        Output folder name for saving cosine similarity results and plots
        
    Returns:
    --------
    pd.DataFrame
        DataFrame containing gene names and cosine similarity scores
    """

    if not os.path.exists(folder):
        os.makedirs(folder)
        print(f"Folder '{folder}' was created.")
    else:
        print(f"Folder '{folder}' already exists.")
    
    # Check if target gene exists
    if target_gene not in adata.var_names:
        raise ValueError(f"Gene '{target_gene}' not found in the dataset")
    
    # Get gene expression matrix
    if issparse(adata.X):
        gene_expression = adata.X.toarray().T  # Transpose: rows = genes, columns = cells
    else:
        gene_expression = adata.X.T
    
    # Get target gene index
    target_idx = np.where(adata.var_names == target_gene)[0][0]
    
    # Calculate cosine similarity
    # Use sklearn's cosine_similarity function
    similarities = cosine_similarity(gene_expression[target_idx:target_idx+1], gene_expression)[0]
    
    # Create results DataFrame
    results = pd.DataFrame({
        'gene': adata.var_names,
        'cosine_similarity': similarities
    })
    
    # Remove target gene itself (similarity = 1)
    results = results[results['gene'] != target_gene]
    
    # Sort by similarity in descending order
    results = results.sort_values('cosine_similarity', ascending=False)
    
    # If top_n is specified, return only top n genes
    if top_n is not None:
        results = results.head(top_n)
    
    # Save results
    results.to_csv(f"{folder}/cosine_similarity_{target_gene}.csv", index=False)
    
    # Plot similarity distribution
    plt.figure(figsize=(10, 6))
    plt.hist(results['cosine_similarity'], bins=50, alpha=0.7, edgecolor='black')
    plt.xlabel('Cosine Similarity')
    plt.ylabel('Number of Genes')
    plt.title(f'Distribution of Cosine Similarity with {target_gene}')
    plt.axvline(results['cosine_similarity'].mean(), color='red', linestyle='--', 
                label=f'Mean: {results["cosine_similarity"].mean():.3f}')
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{folder}/cosine_similarity_distribution_{target_gene}.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # Plot bar chart for top genes
    if top_n is not None and top_n <= 20:  # Only plot bar chart for a small number of genes
        plt.figure(figsize=(12, 8))
        top_results = results.head(top_n)
        plt.barh(range(len(top_results)), top_results['cosine_similarity'])
        plt.yticks(range(len(top_results)), top_results['gene'])
        plt.xlabel('Cosine Similarity')
        plt.title(f'Top {top_n} Genes Most Similar to {target_gene}')
        plt.gca().invert_yaxis()
        plt.tight_layout()
        plt.savefig(f"{folder}/top_{top_n}_similar_genes_{target_gene}.png", dpi=300, bbox_inches='tight')
        plt.close()
    
    print(f"Cosine similarity analysis completed for {target_gene}")
    print(f"Results saved to {folder}/")
    
    return results



def js_divergence(p: np.ndarray, q: np.ndarray) -> float:
    """
    Calculate Jensen-Shannon Divergence (symmetric).
    
    JS(P||Q) = 0.5 * KL(P || M) + 0.5 * KL(Q || M)
    where M = 0.5 * (P + Q)
    
    Parameters
    ----------
    p, q : np.ndarray
        Probability distributions (must sum to 1)
        
    Returns
    -------
    float
        Jensen-Shannon Divergence value (0 to 1, after normalization)
    """
    m = 0.5 * (p + q)
    js = 0.5 * np.sum(kl_div(p, m)) + 0.5 * np.sum(kl_div(q, m))
    # Normalize to [0, 1] range (divide by ln(2) for natural log)
    return float(js / np.log(2))


def _permute_expression(
    gene_expr: np.ndarray,
    coords: np.ndarray,
    method: str,
    *,
    block_n_bins: int = 10,
    rng: np.random.Generator | None = None,
) -> np.ndarray:
    """
    Generate a permuted expression vector under different null models.

    Parameters
    ----------
    gene_expr : np.ndarray
        Expression vector of shape (n_spots,).
    coords : np.ndarray
        Spatial coordinates of shape (n_spots, 2).
    method : str
        Permutation method:
        - 'random': fully random permutation of expression values across spots.
        - 'spatial_blocks': block permutation that preserves within-block local structure.
    block_n_bins : int, default 10
        Number of blocks per axis for 'spatial_blocks' (total blocks = block_n_bins^2).
    rng : np.random.Generator | None
        Random number generator for reproducibility.

    Returns
    -------
    np.ndarray
        Permuted expression vector of shape (n_spots,).
    """
    if rng is None:
        rng = np.random.default_rng()

    method = str(method).lower()
    if method == "random":
        return rng.permutation(gene_expr)

    if method != "spatial_blocks":
        raise ValueError(f"Unknown permutation method: {method}. Use 'random' or 'spatial_blocks'.")

    # Spatial block permutation:
    # 1) assign each spot to a spatial block (tile) on a coarse grid
    # 2) permute the block labels (shuffle blocks)
    # 3) for each original block, fill its spots with values taken from the mapped block
    x = coords[:, 0]
    y = coords[:, 1]

    # Build block grid edges
    x_edges = np.linspace(x.min(), x.max(), block_n_bins + 1)
    y_edges = np.linspace(y.min(), y.max(), block_n_bins + 1)

    # Assign blocks (0..block_n_bins-1 for each axis)
    bx = np.clip(np.digitize(x, x_edges) - 1, 0, block_n_bins - 1)
    by = np.clip(np.digitize(y, y_edges) - 1, 0, block_n_bins - 1)
    block_id = bx * block_n_bins + by  # 0..block_n_bins^2-1

    unique_blocks = np.unique(block_id)
    permuted_blocks = rng.permutation(unique_blocks)
    block_map = dict(zip(unique_blocks.tolist(), permuted_blocks.tolist()))

    permuted_expr = np.empty_like(gene_expr)
    for old_b in unique_blocks:
        new_b = block_map[int(old_b)]
        old_idx = np.where(block_id == old_b)[0]
        new_idx = np.where(block_id == new_b)[0]

        # If mapped block has different number of spots, sample with replacement.
        if len(new_idx) == 0:
            # Fallback: if a block is empty (rare with linspace+digitize), use random permutation
            permuted_expr[old_idx] = rng.permutation(gene_expr[old_idx])
        elif len(new_idx) >= len(old_idx):
            chosen = rng.choice(new_idx, size=len(old_idx), replace=False)
            permuted_expr[old_idx] = gene_expr[chosen]
        else:
            chosen = rng.choice(new_idx, size=len(old_idx), replace=True)
            permuted_expr[old_idx] = gene_expr[chosen]

    return permuted_expr


def filter_genes_by_js_divergence(
    adata,
    marker_gene: str,
    gene_set: list,
    num_permutations: int = 1000,
    pseudocount: float = 1e-10,
    folder: str | None = None,
    fdr_method: str = 'fdr_bh',
    spatial_key: str = "spatial",
    n_bins: int = 20,
    permutation_method: str = "random",
    block_n_bins: int = 10,
    random_state: int | None = None,
) -> pd.DataFrame:
    """
    Filter genes with expression patterns similar to a marker gene using Jensen-Shannon Divergence
    and permutation testing, with multiple testing correction.

    Parameters
    ----------
    adata : AnnData
        Spatial transcriptomics data object
    marker_gene : str
        Marker gene name used as reference.
    gene_set : list
        Gene set to test (automatically removes genes not in the matrix and the marker gene itself).
    num_permutations : int, default 1000
        Number of permutations for testing.
    pseudocount : float, default 1e-10
        Pseudocount to prevent log(0).
    folder : str | None, default None
        If provided, save results to this directory.
    fdr_method : str, default 'fdr_bh'
        Multiple testing correction method, passed to statsmodels.multipletests.
    spatial_key : str, default "spatial"
        Key in adata.obsm containing spatial coordinates (n_spots x 2).
    n_bins : int, default 20
        Number of bins per axis for 2D histogram (higher -> finer grid).
    permutation_method : str, default "random"
        Null model for permutation test:
        - "random": fully random permutation across spots (current behavior)
        - "spatial_blocks": spatial block permutation (preserves within-block local structure)
    block_n_bins : int, default 10
        For permutation_method="spatial_blocks": number of blocks per axis.
    random_state : int | None, default None
        Random seed for reproducibility.

    Returns
    -------
    pd.DataFrame
        Columns: ['gene', 'js_divergence', 'p_value', 'fdr'], sorted by fdr and js_divergence in ascending order.
    """
    if spatial_key not in adata.obsm:
        raise ValueError(f"Spatial coordinates not found in adata.obsm['{spatial_key}']")

    expression_matrix = adata.to_df().T  # genes x cells

    if marker_gene not in expression_matrix.index:
        raise ValueError("Marker gene not found in the expression matrix.")

    # Only keep genes that exist in the gene set and are not the marker gene
    genes_to_test = [g for g in gene_set if g in expression_matrix.index and g != marker_gene]
    if len(genes_to_test) == 0:
        return pd.DataFrame(columns=["gene", "js_divergence", "p_value", "fdr"]) 

    coords = adata.obsm[spatial_key]
    x, y = coords[:, 0], coords[:, 1]
    x_bins = np.linspace(x.min(), x.max(), n_bins + 1)
    y_bins = np.linspace(y.min(), y.max(), n_bins + 1)

    def gene_2d_prob(gene_expr: np.ndarray) -> np.ndarray:
        H, _, _ = np.histogram2d(
            x, y, bins=[x_bins, y_bins], weights=gene_expr
        )
        H = H + pseudocount
        H = H / H.sum()
        return H.ravel()

    q_marker_dist = gene_2d_prob(expression_matrix.loc[marker_gene].values)

    results: list[dict] = []
    rng = np.random.default_rng(random_state)

    for gene in genes_to_test:
        p_dist = gene_2d_prob(expression_matrix.loc[gene].values)
        # Calculate Jensen-Shannon Divergence (symmetric)
        observed_js = js_divergence(p_dist, q_marker_dist)

        # Build null distribution through permutation (based on original expression vector to preserve marginal distribution characteristics)
        null_js_values = []
        original_vec = expression_matrix.loc[gene].values
        for _ in range(num_permutations):
            permuted = _permute_expression(
                original_vec,
                coords,
                permutation_method,
                block_n_bins=block_n_bins,
                rng=rng,
            )
            perm_p = gene_2d_prob(permuted)
            null_js = js_divergence(perm_p, q_marker_dist)
            null_js_values.append(null_js)

        null_js_values = np.asarray(null_js_values)
        extreme_count = int(np.sum(null_js_values <= observed_js))
        p_value = (extreme_count + 1) / (num_permutations + 1)

        results.append({
            "gene": gene,
            "js_divergence": observed_js,
            "p_value": p_value,
        })

    if len(results) == 0:
        return pd.DataFrame(columns=["gene", "js_divergence", "p_value", "fdr"]) 

    results_df = pd.DataFrame(results)
    reject, fdr_values, _, _ = multipletests(results_df["p_value"].values, alpha=0.05, method=fdr_method)
    results_df["fdr"] = fdr_values

    results_df_sorted = results_df.sort_values(by=["fdr", "js_divergence"], ascending=True).reset_index(drop=True)

    if folder is not None:
        if not os.path.exists(folder):
            os.makedirs(folder)
        results_df_sorted.to_csv(os.path.join(folder, f"js_filter_{marker_gene}.csv"), index=False)

    return results_df_sorted



def gene_sets_js_divergence_from_adata(
    adata,
    genesets: dict,
    num_permutations: int = 1000,
    pseudocount: float = 1e-10,
    folder: str = "05.KL",
    fdr_method: str = 'fdr_bh',
    spatial_key: str = "spatial",
    n_bins: int = 20,
    permutation_method: str = "random",
    block_n_bins: int = 10,
    random_state: int | None = None,
) -> pd.DataFrame:
    """
    Based on AnnData and given genesets (dict: celltype -> [genes]),
    for each celltype:
      1) Select the gene with highest average expression in adata from the gene list as marker;
      2) Calculate Jensen-Shannon Divergence of other genes in the list relative to this marker;
      3) Obtain p-value through permutation testing and perform multiple correction to get q-value (fdr).

    Returns a combined results table with columns:
      ['celltype', 'marker', 'gene', 'js_divergence', 'p_value', 'fdr'].
    Also writes CSV files for each celltype and a summary CSV in the folder.
    """
    if not os.path.exists(folder):
        os.makedirs(folder)

    all_rows: list[pd.DataFrame] = []

    for celltype, gene_list in genesets.items():
        # Filter to keep only genes present in the data
        available_genes = [g for g in gene_list if g in adata.var_names]
        if len(available_genes) < 2:
            # Need at least 2 genes (one as marker, others for comparison)
            continue

        # Select marker: from the gene list, choose the one with highest average expression (across all cells/samples)
        exp_df = adata.to_df().T  # genes x cells
        sub_mat = exp_df.loc[available_genes]
        mean_expr = sub_mat.mean(axis=1)
        marker_gene = str(mean_expr.idxmax())

        # Compare other genes with marker using JS divergence
        compare_genes = [g for g in available_genes if g != marker_gene]

        # Use existing JS divergence function for statistical testing
        res = filter_genes_by_js_divergence(
            adata=adata,
            marker_gene=marker_gene,
            gene_set=compare_genes,
            num_permutations=num_permutations,
            pseudocount=pseudocount,
            folder=None,
            fdr_method=fdr_method,
            spatial_key=spatial_key,
            n_bins=n_bins,
            permutation_method=permutation_method,
            block_n_bins=block_n_bins,
            random_state=random_state,
        )

        if res.empty:
            continue

        res.insert(0, 'marker', marker_gene)
        res.insert(0, 'celltype', celltype)

        # Save results for each celltype
        res.to_csv(os.path.join(folder, f"jsd_{celltype}.csv"), index=False)
        all_rows.append(res)

    if len(all_rows) == 0:
        return pd.DataFrame(columns=['celltype', 'marker', 'gene', 'js_divergence', 'p_value', 'fdr'])

    combined = pd.concat(all_rows, axis=0, ignore_index=True)
    combined.to_csv(os.path.join(folder, "jsd_all_celltypes.csv"), index=False)
    return combined
