"""
STEA: Spatial Transcriptome Enrichment Analysis
A Python package for deconvolution of spatial transcriptomics datasets.
"""

import os
from pickle import TRUE
from .utils import (
    qcplot,
    marker_rank_analysis,
    stgsea,
    calculate_gene_cosine_similarity,
    filter_genes_by_js_divergence,
    gene_sets_js_divergence_from_adata,
)

__version__ = "0.1.0"
__author__ = "LI QIAN"
__email__ = "liqian0901@connect.hku.hk"

# Main workflow function
def stea(adata, genesets, QC = TRUE, folder="0.STEA"):
    """
    STEA analysis pipeline
    
    Parameters:
    -----------
    adata : AnnData
    genesets : dict, formats {celltype: [marker list]}
    QC : bool, default=True, whether to perform quality control analysis
    folder : str, default="temp", output folder name
        
    Returns:
    --------
    dict
        Dictionary containing analysis results
    """
    if not os.path.exists(folder):
        os.makedirs(folder)
        print(f"Folder '{folder}' was created.")
    else:
        print(f"Folder '{folder}' already exists.")
     
    # Save current working directory
    original_dir = os.getcwd()
    try:
        # Switch to output directory
        os.chdir(folder)
        
        print("="*60)
        print("STEA Analysis Pipeline Started")
        print("="*60)
        
        # Run QC inside the STEA output directory
        if QC:
            print("\n1. Quality Control Analysis...")
            qcplot(adata, folder="01.QC")
        else:
            print("Quality Control Analysis skipped.")
        
        # 2. Marker gene ranking analysis
        print("\n1. Marker Gene Ranking Analysis...")
        marker_results = marker_rank_analysis(adata, genesets, folder="02.MarkerRank")
        
        # 3. Gene set enrichment analysis
        print("\n2. Gene Set Enrichment Analysis...")
        es_results = stgsea(adata, genesets, folder="03.ES")
        
        print("\n" + "="*60)
        print("STEA Analysis Pipeline Completed Successfully!")
        print("="*60)
        
        return  es_results
        
    finally:
        # Restore original working directory
        os.chdir(original_dir)


def cosine_similarity(adata, gene, top_n=None, folder="04.CosineSimilarity"):
    """
    Calculate cosine similarity between a gene and all other genes in the dataset

    Parameters:
    -----------
    adata : AnnData
    gene : str
    top_n : int, optional, default=None
    folder : str, default="04.CosineSimilarity"
        
    """
    return calculate_gene_cosine_similarity(adata, gene, top_n, folder)

def js_divergence(
    adata,
    genesets,
    num_permutations=1000,
    pseudocount=1e-10,
    folder="05.JS",
    fdr_method='fdr_bh',
    spatial_key="spatial",
    n_bins=20,
    permutation_method="random",
    block_n_bins=10,
    random_state=None,
):
    """
    Calculate Jensen-Shannon Divergence for gene sets using AnnData object.
    
    Note: This function uses Jensen-Shannon Divergence (symmetric) for
    spatial expression pattern comparison.

    Parameters:
    -----------
    adata : AnnData
        Spatial transcriptomics data object
    genesets : dict
        Gene sets dictionary with format {celltype: [marker list]}
    num_permutations : int, default=1000
        Number of permutations for testing
    pseudocount : float, default=1e-10
        Pseudocount to prevent log(0)
        folder : str, default="05.JS"
        Output folder name
    fdr_method : str, default='fdr_bh'
        Multiple testing correction method
    spatial_key : str, default "spatial"
        Key in adata.obsm containing spatial coordinates
    n_bins : int, default 20
        Number of bins per axis for 2D histogram
    permutation_method : str, default "random"
        Null model for permutation test: "random" or "spatial_blocks"
    block_n_bins : int, default 10
        For permutation_method="spatial_blocks": number of blocks per axis
    random_state : int | None, default None
        Random seed for reproducibility
        
    Returns:
    --------
    pd.DataFrame
        Combined results table with columns: ['celltype', 'marker', 'gene', 'js_divergence', 'p_value', 'fdr']
        Note: 'js_divergence' contains Jensen-Shannon Divergence values (0 to 1, where 0 = identical, 1 = maximally different)
    """
    return gene_sets_js_divergence_from_adata(
        adata,
        genesets,
        num_permutations,
        pseudocount,
        folder,
        fdr_method,
        spatial_key,
        n_bins,
        permutation_method,
        block_n_bins,
        random_state,
    )


# Export main functions and classes
__all__ = [
    'stea',
    'qcplot',
    'marker_rank_analysis',
    'stgsea',
    'calculate_gene_cosine_similarity',
    'cosine_similarity',
    'filter_genes_by_js_divergence',
    'js_divergence',
    'gene_sets_js_divergence_from_adata',
]

