#!/usr/bin/env python3
"""
Configuration parameters for read clustering and polishing.

This module contains all tunable parameters for the read clustering pipeline,
including BLAST parameters, VSEARCH parameters, clustering settings, and
polishing parameters.
"""

# =============================================================================
# BLAST PARAMETERS - NOT NEEDED IF VSEARCH
# =============================================================================

BLAST_PARAMS = {
    'identity': 80,  # % identity for clustering. Default 90%
    'query_coverage': 80,  # 95% query coverage threshold. Default 90%. 
    'evalue': 1e-10,  # Strict E-value for high-quality matches
    'word_size': 7,  # Standard word size
    'dust': 'no',  # Disable dust filtering
    'soft_masking': 'false',
    'gapopen': 5,  # Gap opening penalty
    'gapextend': 2,  # Gap extension penalty
    'penalty': -1,  # Mismatch penalty
    'reward': 1,  # Match reward
    'max_target_seqs': 1000000,  # Allow multiple hits per query. Default 100000
    'task': 'blastn'  # Standard BLASTN
}

# =============================================================================
# BLAST CLUSTERING PARAMETERS - NOT NEEDED IF VSEARCH
# =============================================================================

CLUSTERING_PARAMS = {
    'coverage_threshold': 95,  # Minimum coverage percentage to include in cluster. Default 90
    'identity_threshold': 90,  # Minimum identity percentage for clustering. Default 90
    'min_cluster_size': 5,  # Minimum number of reads in a cluster. Default 5
    'max_clusters': None,  # Maximum number of clusters (None = unlimited)
    'max_read_length': None,  # Maximum read length to consider (None = unlimited)
    'use_fast_clustering': False,  # Use VSEARCH instead of BLAST for clustering
}

# =============================================================================
# VSEARCH CLUSTERING PARAMETERS
# =============================================================================

VSEARCH_PARAMS = {
    # Clustering algorithm options:
    # - cluster_fast: Fastest, good for similar-length sequences
    # - cluster_smallmem: Better quality, moderate memory usage
    # - cluster_size: Best quality, handles diverse lengths well
    # - cluster_unoise: Denoising algorithm, highest quality but slowest
    'clustering_algorithm': 'cluster_size',
    
    # Identity threshold (0.80-0.95, higher = stricter clustering)
    # - 0.80: Very inclusive, many sequences per cluster
    # - 0.85: Inclusive, good for noisy data
    # - 0.90: Balanced, recommended for most cases
    # - 0.95: Strict, high-quality clusters
    'identity_threshold': 0.9, #0.9 is good
    
    # Coverage thresholds (0.70-0.90, higher = stricter)
    # - query_coverage: How much of the query sequence must align
    # - target_coverage: How much of the target sequence must align
    # - 0.80 is recommended for good quality without being too strict
    'query_coverage': 0.9, #0.95 was good , qnrs15 used 0.9*
    'target_coverage': 0.85, #0.95 was good, qnrs15 used 0.9*
    
    # Strand handling: 'both', 'plus', or 'minus'
    'strand': 'both',
    
    # Size annotations for better clustering
    'sizein': True,   # Input contains size annotations
    'sizeout': True,  # Output contains size annotations
}

# =============================================================================
# MAFFT CONSENSUS POLISHING PARAMETERS - NO TNEEDED IF RACON
# =============================================================================

POLISHING_PARAMS = {
    'min_coverage_threshold': 0.3,  # Minimum coverage for consensus (50% of reads must have a base). Default 0.5
    'gap_threshold': 0.3,  # Threshold for calling gaps (1.0 = 100% of reads must have gap, 0.5 = 50% of reads must have gap) %0.1 is good
    'ambiguous_threshold': 1.0,  # Threshold for calling ambiguous bases (1.0 = always call most frequent base, never N) %1 is good
    'max_reads_for_polishing': 50,  # Maximum number of reads to use for polishing large clusters %50 is good
    # MAFFT alignment parameters optimized for 95% identical reads with varying sizes (trimmed ends)
    'mafft_maxiterate': 1000,  # Maximum iterations for MAFFT (high for optimal alignment of very similar sequences)
    'mafft_retree': 2,  # Number of tree rebuilding cycles (better alignment quality)
    'mafft_op': 1.53,  # Gap opening penalty (moderate for similar sequences)
    'mafft_ep': 0.123,  # Gap extension penalty (low for small indels and length differences)
    'mafft_use_localpair': False,  # Use global alignment for full sequence coverage
    'mafft_threads': 1,  # Number of threads (single thread for consistency)
}

# =============================================================================
# RACON POLISHING PARAMETERS
# =============================================================================

RACON_PARAMS = {
    # RACON consensus parameters
    'match': 3,  # Match score (default: 3)
    'mismatch': -5,  # Mismatch penalty (default: -5)
    'gap': -4,  # Gap penalty (default: -4)
    'window_length': 500,  # Window length for consensus (default: 500)
    'quality_threshold': 10.0,  # Quality threshold for consensus (default: 10.0)
    'error_threshold': 0.3,  # Error threshold for consensus (default: 0.3)
    'threads': 1,  # Number of threads (default: 1)
    'include_unpolished': False,  # Include unpolished regions (default: False) *??
    'split': False,  # Split reads at gaps (default: False)
    'fragments': False,  # Use fragments (default: False)
    'max_reads_for_polishing': 50,  # Maximum number of reads to use for polishing large clusters
}

# =============================================================================
# REDUNDANT CONSENSUS REMOVAL PARAMETERS
# =============================================================================

REDUNDANT_REMOVAL_PARAMS = {
    'enabled': False,  # Enable redundant consensus sequence removal
    'identity_threshold': 0.95,  # Identity threshold for considering sequences redundant (0.95 = 95%)
    'query_coverage': 0.80,  # Query coverage threshold for redundancy detection
    'target_coverage': 0.80,  # Target coverage threshold for redundancy detection
    'strand': 'both',  # Strand handling: 'both', 'plus', or 'minus'
    'output_file_suffix': '_deduplicated',  # Suffix for deduplicated output file
}

# =============================================================================
# DEBUG AND OUTPUT PARAMETERS
# =============================================================================

DEBUG_PARAMS = {
    'save_blast_outputs': False,  # Save BLAST output files for debugging
    'verbose': False,  # Print detailed output for each operation
    'save_intermediate': False,  # Save intermediate files
    'output_longest_reads': True,  # Output longest read from each cluster to separate FASTA file
}

# =============================================================================
# ALGORITHM SELECTION
# =============================================================================

# Available clustering algorithms
CLUSTERING_ALGORITHMS = {
    'blast': 'cluster_by_blast',
    'vsearch': 'cluster_by_vsearch'
}

# Available polishing algorithms
POLISHING_ALGORITHMS = {
    'mafft': 'polish_by_mafft',
    'racon': 'polish_by_racon',
    'simple': 'polish_simple'  # Simple consensus without alignment
}

# Default algorithm selections
DEFAULT_CLUSTERING_ALGORITHM = 'blast'
DEFAULT_POLISHING_ALGORITHM = 'mafft'
