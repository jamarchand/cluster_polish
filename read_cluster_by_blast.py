#!/usr/bin/env python3
"""
Read Clustering by BLAST Tool

This script implements a greedy clustering algorithm that groups reads based on 
high sequence similarity using BLAST. It picks the longest read as a representative,
finds all reads that cover >95% of it, and creates clusters.

Author: Assistant
Date: 2024
"""

import os
import sys
import subprocess
import tempfile
import shutil
import random
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import AlignIO
import argparse
import logging
import re
from collections import Counter

# =============================================================================
# TUNABLE PARAMETERS - Modify these as needed
# =============================================================================

# BLAST parameters for read clustering
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

# Clustering parameters
CLUSTERING_PARAMS = {
    'coverage_threshold': 95,  # Minimum coverage percentage to include in cluster. Default 90
    'identity_threshold': 90,  # Minimum identity percentage for clustering. Default 90
    'min_cluster_size': 5,  # Minimum number of reads in a cluster. Default 5
    'max_clusters': None,  # Maximum number of clusters (None = unlimited)
    'max_read_length': None,  # Maximum read length to consider (None = unlimited)
    'use_fast_clustering': False,  # Use VSEARCH instead of BLAST for clustering
}

# VSEARCH clustering parameters (when use_fast_clustering=True)
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
    'identity_threshold': 0.90, #0.9 is good
    
    # Coverage thresholds (0.70-0.90, higher = stricter)
    # - query_coverage: How much of the query sequence must align
    # - target_coverage: How much of the target sequence must align
    # - 0.80 is recommended for good quality without being too strict
    'query_coverage': 0.90, #0.95 was good , qnrs15 used 0.9*
    'target_coverage': 0.90, #0.95 was good, qnrs15 used 0.9*
    
    # Strand handling: 'both', 'plus', or 'minus'
    'strand': 'both',
    
    # Size annotations for better clustering
    'sizein': True,   # Input contains size annotations
    'sizeout': True,  # Output contains size annotations
}

# Consensus polishing parameters
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
#best 14 
# Debug and output parameters
DEBUG_PARAMS = {
    'save_blast_outputs': False,  # Save BLAST output files for debugging
    'verbose': False,  # Print detailed output for each operation
    'save_intermediate': False,  # Save intermediate files
}

# =============================================================================
# END TUNABLE PARAMETERS
# =============================================================================

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class ReadClusterer:
    def __init__(self, input_file, output_dir, coverage_threshold=None, identity_threshold=None,
                 min_cluster_size=None, max_clusters=None, max_read_length=None, save_blast_outputs=None, 
                 verbose=None, save_intermediate=None, polish_clusters=None, use_fast_clustering=None):
        """
        Initialize the read clusterer with input file and parameters.
        
        Args:
            input_file (str): Path to input FASTA file
            output_dir (str): Directory for output files
            coverage_threshold (int): Minimum coverage percentage for clustering
            identity_threshold (int): Minimum identity percentage for clustering
            min_cluster_size (int): Minimum reads in a cluster
            max_clusters (int): Maximum number of clusters
            save_blast_outputs (bool): Save BLAST output files for debugging
            verbose (bool): Print detailed output
            save_intermediate (bool): Save intermediate files
            polish_clusters (bool): Whether to calculate consensus sequences for clusters
            use_fast_clustering (bool): Whether to use VSEARCH instead of BLAST for clustering
        """
        self.input_file = Path(input_file)
        self.output_dir = Path(output_dir)
        
        # Use provided parameters or defaults
        self.coverage_threshold = coverage_threshold if coverage_threshold is not None else CLUSTERING_PARAMS['coverage_threshold']
        self.identity_threshold = identity_threshold if identity_threshold is not None else CLUSTERING_PARAMS['identity_threshold']
        self.min_cluster_size = min_cluster_size if min_cluster_size is not None else CLUSTERING_PARAMS['min_cluster_size']
        self.max_clusters = max_clusters if max_clusters is not None else CLUSTERING_PARAMS['max_clusters']
        self.max_read_length = max_read_length if max_read_length is not None else CLUSTERING_PARAMS['max_read_length']
        self.save_blast_outputs = save_blast_outputs if save_blast_outputs is not None else DEBUG_PARAMS['save_blast_outputs']
        self.verbose = verbose if verbose is not None else DEBUG_PARAMS['verbose']
        self.save_intermediate = save_intermediate if save_intermediate is not None else DEBUG_PARAMS['save_intermediate']
        self.polish_clusters = polish_clusters if polish_clusters is not None else False
        self.use_fast_clustering = use_fast_clustering if use_fast_clustering is not None else CLUSTERING_PARAMS['use_fast_clustering']
        
        # Use parameters from the top of the script
        self.blast_params = BLAST_PARAMS.copy()
        self.polishing_params = POLISHING_PARAMS.copy()
        
        # Validate input file and create output directory
        self._validate_inputs()
        self.output_dir.mkdir(exist_ok=True)
        
        # Initialize clustering statistics
        self.cluster_stats = {
            'total_reads': 0,
            'total_clusters': 0,
            'clusters_created': 0,
            'reads_clustered': 0,
            'reads_unclustered': 0,
            'clusters_polished': 0
        }
        
        # Initialize cluster list for tracking parent sequences
        self.cluster_list = []
        
        # Initialize growing cluster list file
        self.growing_cluster_file = self.output_dir / "cluster_list_growing.fasta"
        
        # Initialize polished consensus file
        if self.polish_clusters:
            self.consensus_file = self.output_dir / "cluster_consensus.fasta"
            self.consensus_list = []
        
    def _validate_inputs(self):
        """Validate that input file exists and is readable."""
        if not self.input_file.exists():
            raise FileNotFoundError(f"Input file not found: {self.input_file}")
        
        # Check if input is FASTA
        self.input_format = self._detect_format(self.input_file)
        if self.input_format != 'fasta':
            raise ValueError(f"Input file must be in FASTA format, detected: {self.input_format}")
        
        # Note: MAFFT is required for consensus calculation
        if self.polish_clusters:
            logger.info("Nanopore consensus calculation enabled (MAFFT required for alignment)")
        
        logger.info(f"Input format: {self.input_format}")
        
    def _detect_format(self, file_path):
        """Detect if file is FASTA or FASTQ format."""
        with open(file_path, 'r') as f:
            first_line = f.readline().strip()
            if first_line.startswith('@'):
                return 'fastq'
            elif first_line.startswith('>'):
                return 'fasta'
            else:
                raise ValueError(f"Unknown file format: {file_path}")
    
    def _sort_reads_by_length(self, input_file):
        """Sort reads by length (longest first) and return as list."""
        reads = list(SeqIO.parse(input_file, 'fasta'))
        
        # Filter by max read length if specified
        if self.max_read_length is not None:
            original_count = len(reads)
            reads = [read for read in reads if len(read.seq) <= self.max_read_length]
            filtered_count = len(reads)
            logger.info(f"Filtered reads by max length {self.max_read_length}: {filtered_count}/{original_count} reads kept")
        
        reads.sort(key=lambda x: len(x.seq), reverse=True)
        
        self.cluster_stats['total_reads'] = len(reads)
        logger.info(f"Sorted {len(reads)} reads by length (longest first)")
        
        if self.verbose:
            for i, read in enumerate(reads[:5]):  # Show top 5 longest reads
                logger.debug(f"  {i+1}. {read.id}: {len(read.seq)} bp")
        
        return reads
    
    def _fast_clustering_vsearch(self, reads, temp_dir):
        """Use VSEARCH for fast clustering instead of BLAST."""
        try:
            # Check if VSEARCH is available
            result = subprocess.run(['vsearch', '--version'], capture_output=True, text=True, check=True)
            logger.info(f"VSEARCH version: {result.stdout.strip()}")
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            logger.error(f"VSEARCH is required for fast clustering but not found: {e}")
            logger.error("Please install VSEARCH: conda install -c bioconda vsearch")
            raise RuntimeError("VSEARCH is required for fast clustering. Please install it first.")
        
        logger.info("Using VSEARCH for fast clustering")
        
        # Create temporary input file
        input_file = temp_dir / "reads_for_vsearch.fasta"
        SeqIO.write(reads, input_file, 'fasta')
        
        # VSEARCH clustering output files
        output_prefix = temp_dir / "vsearch_clusters"
        uc_file = temp_dir / "clusters.uc"
        
        # Run VSEARCH clustering with configurable parameters for better quality
        # Algorithm and thresholds are configurable via VSEARCH_PARAMS at top of script
        cmd = [
            'vsearch', f'--{VSEARCH_PARAMS["clustering_algorithm"]}', str(input_file),
            '--id', str(VSEARCH_PARAMS['identity_threshold']),
            '--strand', VSEARCH_PARAMS['strand'],
            '--uc', str(uc_file),
            '--centroids', str(output_prefix) + '.fasta'
        ]
        
        # Add coverage thresholds if specified
        if VSEARCH_PARAMS['query_coverage']:
            cmd.extend(['--query_cov', str(VSEARCH_PARAMS['query_coverage'])])
        if VSEARCH_PARAMS['target_coverage']:
            cmd.extend(['--target_cov', str(VSEARCH_PARAMS['target_coverage'])])
        
        # Add size annotations if enabled
        if VSEARCH_PARAMS['sizein']:
            cmd.append('--sizein')
        if VSEARCH_PARAMS['sizeout']:
            cmd.append('--sizeout')
        
        logger.info(f"Running VSEARCH clustering: {' '.join(cmd)}")
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        logger.info("VSEARCH clustering completed")
        
        # Parse VSEARCH output
        clusters = self._parse_vsearch_output(uc_file, reads)
        
        return clusters
    
    def _parse_vsearch_output(self, uc_file, reads):
        """Parse VSEARCH UC output format and return clusters."""
        clusters = []
        cluster_id_to_index = {}  # Map VSEARCH cluster ID to cluster index
        
        # Read UC file
        with open(uc_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                fields = line.split('\t')
                if len(fields) < 10:
                    continue
                
                record_type = fields[0]
                cluster_id = fields[1]  # VSEARCH cluster identifier
                read_id = fields[8]     # Read identifier
                
                if record_type == 'S':  # Seed sequence (representative)
                    # Create new cluster
                    read_obj = next((r for r in reads if r.id == read_id), None)
                    if read_obj:
                        clusters.append([read_obj])
                        cluster_id_to_index[cluster_id] = len(clusters) - 1  # Map cluster ID to cluster index
                
                elif record_type == 'H':  # Hit sequence (member)
                    # Add to existing cluster
                    read_obj = next((r for r in reads if r.id == read_id), None)
                    if read_obj and cluster_id in cluster_id_to_index:  # Use cluster ID to find cluster
                        cluster_idx = cluster_id_to_index[cluster_id]
                        clusters[cluster_idx].append(read_obj)
        
        # Filter clusters by minimum size and quality
        filtered_clusters = []
        for cluster in clusters:
            if len(cluster) >= self.min_cluster_size:
                # For VSEARCH clusters, filter out sequences that are too divergent
                if self.use_fast_clustering:
                    filtered_cluster = self._filter_vsearch_cluster_quality(cluster)
                    if len(filtered_cluster) >= self.min_cluster_size:
                        filtered_clusters.append(filtered_cluster)
                else:
                    filtered_clusters.append(cluster)
        
        # Debug: Show cluster size distribution
        if clusters:
            cluster_sizes = [len(cluster) for cluster in clusters]
            filtered_sizes = [len(cluster) for cluster in filtered_clusters]
            logger.info(f"Original cluster size distribution: min={min(cluster_sizes)}, max={max(cluster_sizes)}, mean={sum(cluster_sizes)/len(cluster_sizes):.1f}")
            if filtered_sizes:
                logger.info(f"Filtered cluster size distribution: min={min(filtered_sizes)}, max={max(filtered_sizes)}, mean={sum(filtered_sizes)/len(filtered_sizes):.1f}")
            logger.info(f"Clusters with size >= {self.min_cluster_size}: {len(filtered_clusters)}/{len(clusters)}")
            
            # Show some examples of small clusters
            small_clusters = [cluster for cluster in clusters if len(cluster) < self.min_cluster_size]
            if small_clusters:
                logger.info(f"Examples of small clusters (size < {self.min_cluster_size}):")
                for i, cluster in enumerate(small_clusters[:5]):  # Show first 5 small clusters
                    logger.info(f"  Small cluster {i+1}: {len(cluster)} reads")
        
        logger.info(f"VSEARCH created {len(clusters)} clusters, {len(filtered_clusters)} meet minimum size and quality requirements")
        
        return filtered_clusters
    
    def _filter_vsearch_cluster_quality(self, cluster):
        """Filter VSEARCH cluster to remove sequences that are too divergent for good MAFFT alignment."""
        if len(cluster) <= 2:
            return cluster  # Too small to filter
        
        # Use the longest sequence as reference
        reference = max(cluster, key=lambda x: len(x.seq))
        reference_length = len(reference.seq)
        
        # Calculate pairwise similarities and filter out divergent sequences
        filtered_cluster = [reference]  # Always keep reference
        
        for read in cluster:
            if read.id == reference.id:
                continue
            
            # Simple length-based similarity check
            length_diff = abs(len(read.seq) - reference_length)
            length_similarity = 1.0 - (length_diff / max(len(read.seq), reference_length))
            
            # Keep sequences that are reasonably similar in length (within 20%)
            if length_similarity >= 0.8:
                filtered_cluster.append(read)
            else:
                logger.debug(f"Filtered out divergent sequence {read.id}: length similarity {length_similarity:.2f}")
        
        logger.debug(f"Cluster quality filtering: {len(cluster)} -> {len(filtered_cluster)} sequences")
        return filtered_cluster
    
    def _create_blast_db(self, query_file, db_name, temp_dir):
        """Create a BLAST database from query sequence."""
        db_base_path = temp_dir / db_name
        
        logger.info(f"Creating BLAST database: {db_base_path}")
        
        cmd = [
            'makeblastdb',
            '-in', str(query_file),
            '-dbtype', 'nucl',
            '-out', str(db_base_path)
        ]
        
        try:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            logger.info(f"Created BLAST database: {db_base_path}")
            return db_base_path
        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to create BLAST database: {e}")
            raise
    
    def _run_blast(self, query_file, db_path, output_file, cluster_id=None):
        """Run BLAST search."""
        cmd = [
            'blastn',
            '-query', str(query_file),
            '-db', str(db_path),
            '-out', str(output_file),
            '-outfmt', '10 qseqid sseqid pident qcovs qstart qend sstart send evalue bitscore sstrand'
        ]
        
        logger.info(f"Running BLAST search: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True)
            logger.info(f"BLAST search completed: {output_file}")
            
            # Save BLAST output to local directory if requested
            if self.save_blast_outputs and cluster_id:
                blast_outputs_dir = self.output_dir / "blast_outputs"
                blast_outputs_dir.mkdir(exist_ok=True)
                
                local_blast_output = blast_outputs_dir / f"blast_cluster_{cluster_id}_output.txt"
                shutil.copy2(output_file, local_blast_output)
                logger.info(f"Saved BLAST output to: {local_blast_output}")
                
        except subprocess.CalledProcessError as e:
            logger.error(f"BLAST search failed: {e}")
            raise
    
    def _parse_blast_results(self, blast_output, query_length):
        """Parse BLAST results and identify reads that meet clustering criteria."""
        hits = {}
        
        if not os.path.exists(blast_output) or os.path.getsize(blast_output) == 0:
            logger.warning(f"BLAST output file is empty or doesn't exist: {blast_output}")
            return hits
        
        logger.info(f"Parsing BLAST results from: {blast_output}")
        
        with open(blast_output, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                
                # Parse CSV format
                fields = line.split(',')
                if len(fields) >= 10:
                    qseqid = fields[0]  # read name
                    sseqid = fields[1]  # query name (should be the representative)
                    pident = float(fields[2])
                    qcovs = float(fields[3])
                    qstart = int(fields[4])  # read start
                    qend = int(fields[5])    # read end
                    sstart = int(fields[6])  # query start
                    send = int(fields[7])    # query end
                    evalue = float(fields[8])
                    sstrand = fields[10] if len(fields) > 10 else 'plus'
                    
                    # Calculate coverage of the query (representative read)
                    query_coverage = abs(send - sstart) + 1
                    query_coverage_percent = (query_coverage / query_length) * 100
                    
                    # Filter by identity and coverage thresholds
                    if pident >= self.identity_threshold and query_coverage_percent >= self.coverage_threshold:
                        if qseqid not in hits:
                            hits[qseqid] = []
                        
                        hits[qseqid].append({
                            'query': sseqid,
                            'start': qstart,
                            'end': qend,
                            'identity': pident,
                            'coverage': query_coverage_percent,
                            'evalue': evalue,
                            'strand': sstrand
                        })
                        
                        if self.verbose:
                            logger.debug(f"Hit for read {qseqid}: identity={pident}%, coverage={query_coverage_percent:.1f}%")
        
        logger.info(f"Found {len(hits)} reads meeting clustering criteria")
        return hits
    
    def _create_cluster(self, representative_read, cluster_reads, cluster_id):
        """Create a cluster file with the representative and cluster reads."""
        cluster_filename = f"cluster_{cluster_id:04d}.fasta"
        cluster_filepath = self.output_dir / cluster_filename
        
        # Write representative read first, then cluster reads
        all_reads = [representative_read] + cluster_reads
        SeqIO.write(all_reads, cluster_filepath, 'fasta')
        
        # Add representative to cluster list for tracking
        representative_copy = SeqRecord(
            seq=representative_read.seq,
            id=f"cluster_{cluster_id:04d}_representative",
            description=f"Representative for cluster {cluster_id} (original: {representative_read.id})"
        )
        self.cluster_list.append(representative_copy)
        
        # Write to growing cluster file immediately
        with open(self.growing_cluster_file, 'a') as f:
            SeqIO.write(representative_copy, f, 'fasta')
        
        logger.info(f"Added cluster {cluster_id} representative to growing list: {self.growing_cluster_file}")
        
        logger.info(f"Created cluster {cluster_id}: {len(all_reads)} reads -> {cluster_filepath}")
        
        return cluster_filepath
    
    def _remove_clustered_reads(self, reads, clustered_read_ids):
        """Remove clustered reads from the reads list."""
        remaining_reads = [read for read in reads if read.id not in clustered_read_ids]
        removed_count = len(reads) - len(remaining_reads)
        
        logger.info(f"Removed {removed_count} clustered reads, {len(remaining_reads)} remaining")
        
        return remaining_reads
    

    
    def _run_mafft_nanopore(self, cluster_file, temp_dir, cluster_id):
        """Run MAFFT with Nanopore-optimized settings for error correction."""
        try:
            # Check if MAFFT is available
            result = subprocess.run(['mafft', '--version'], capture_output=True, text=True, check=True)
            logger.info(f"MAFFT version: {result.stdout.strip()}")
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            logger.warning(f"MAFFT not available, skipping: {e}")
            return None
        
        # Create aligned output file
        aligned_file = temp_dir / f"mafft_aligned_cluster_{cluster_id}.fasta"
        
        # Run MAFFT alignment with Nanopore-optimized settings
        try:
            # Create a completely clean temporary directory with simple names
            clean_temp_dir = Path("/tmp") / f"mafft_cluster_{cluster_id}_{os.getpid()}"
            clean_temp_dir.mkdir(exist_ok=True)
            
            temp_input_file = clean_temp_dir / f"input.fasta"
            temp_output_file = clean_temp_dir / f"output.fasta"
            
            # Read all sequences from cluster file
            reads = list(SeqIO.parse(cluster_file, 'fasta'))
            total_reads = len(reads)
            
            # Randomly sample reads if cluster is too large
            if total_reads > self.polishing_params['max_reads_for_polishing']:
                random.seed(42)  # For reproducible results
                sampled_reads = random.sample(reads, self.polishing_params['max_reads_for_polishing'])
                logger.info(f"Cluster {cluster_id}: Sampling {len(sampled_reads)} reads from {total_reads} total reads for MAFFT")
                reads_to_align = sampled_reads
            else:
                reads_to_align = reads
                logger.info(f"Cluster {cluster_id}: Using all {total_reads} reads for MAFFT")
            
            # Write sampled reads to clean temp directory
            SeqIO.write(reads_to_align, temp_input_file, 'fasta')
            
            # MAFFT settings optimized for 95% identical reads with varying sizes (trimmed ends)
            # Build command with configurable parameters
            cmd = ['mafft']
            
            # Use global or local alignment based on parameter
            if self.polishing_params['mafft_use_localpair']:
                cmd.append('--localpair')
            else:
                cmd.append('--globalpair')
            
            # Add configurable parameters
            cmd.extend(['--maxiterate', str(self.polishing_params['mafft_maxiterate'])])
            cmd.extend(['--retree', str(self.polishing_params['mafft_retree'])])
            cmd.extend(['--op', str(self.polishing_params['mafft_op'])])
            cmd.extend(['--ep', str(self.polishing_params['mafft_ep'])])
            cmd.extend(['--thread', str(self.polishing_params['mafft_threads'])])
            cmd.append('--quiet')
            cmd.append(str(temp_input_file))
            logger.info(f"Running MAFFT with Nanopore settings: {' '.join(cmd)}")
            logger.info(f"Clean input file: {temp_input_file} ({len(reads_to_align)} reads)")
            
            # Check if input file exists and has content
            if not os.path.exists(temp_input_file):
                raise FileNotFoundError(f"Input file does not exist: {temp_input_file}")
            
            if os.path.getsize(temp_input_file) == 0:
                raise ValueError(f"Input file is empty: {temp_input_file}")
            
            # Run MAFFT and capture output
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            logger.info(f"MAFFT alignment completed for cluster {cluster_id}")
            
            # Write MAFFT output to file
            with open(temp_output_file, 'w') as f:
                f.write(result.stdout)
            
            # Copy result back to original temp directory
            shutil.copy2(temp_output_file, aligned_file)
            
            # Clean up clean temp directory
            shutil.rmtree(clean_temp_dir, ignore_errors=True)
            
            return aligned_file
            
        except subprocess.CalledProcessError as e:
            logger.error(f"MAFFT alignment failed for cluster {cluster_id}: {e}")
            if hasattr(e, 'stderr'):
                logger.error(f"MAFFT stderr: {e.stderr}")
            
            # Clean up clean temp directory on error
            if 'clean_temp_dir' in locals():
                shutil.rmtree(clean_temp_dir, ignore_errors=True)
            
            return None
        except Exception as e:
            logger.error(f"Unexpected error during MAFFT alignment for cluster {cluster_id}: {e}")
            
            # Clean up clean temp directory on error
            if 'clean_temp_dir' in locals():
                shutil.rmtree(clean_temp_dir, ignore_errors=True)
            
            return None
    
    def _create_nanopore_consensus(self, cluster_file, cluster_id):
        """Create consensus for Nanopore reads without introducing excessive gaps."""
        try:
            reads = list(SeqIO.parse(cluster_file, 'fasta'))
            if len(reads) < 2:
                logger.warning(f"Cluster {cluster_id} has only {len(reads)} sequence(s), cannot calculate consensus")
                return None
            
            logger.info(f"Creating Nanopore consensus for cluster {cluster_id} with {len(reads)} reads")
            
            # Find the longest sequence length
            max_length = max(len(read.seq) for read in reads)
            logger.info(f"Cluster {cluster_id}: Max length: {max_length}")
            
            # Show sequence lengths for debugging
            for i, read in enumerate(reads):
                logger.debug(f"  Read {i+1}: {len(read.seq)} bp")
            
            # Pad shorter sequences to match the longest (add N's instead of gaps)
            padded_reads = []
            for read in reads:
                if len(read.seq) < max_length:
                    # Pad with N's to match longest sequence
                    padding = 'N' * (max_length - len(read.seq))
                    padded_seq = str(read.seq) + padding
                    padded_reads.append(padded_seq)
                else:
                    padded_reads.append(str(read.seq))
            
            # Calculate consensus position by position using quality-weighted voting
            consensus_seq = ""
            for pos in range(max_length):
                bases = [seq[pos] for seq in padded_reads]
                
                # Count N's and gaps for gap threshold calculation
                n_count = sum(1 for base in bases if base.upper() == 'N')
                gap_count = sum(1 for base in bases if base.upper() in ['-', '.'])
                total_gap_like = n_count + gap_count
                gap_fraction = total_gap_like / len(bases)
                
                # Filter out N's and gaps for base counting
                valid_bases = [base for base in bases if base.upper() not in ['N', '-', '.']]
                
                # Check if gap threshold is met (considering both N's and gaps)
                if gap_fraction >= self.polishing_params['gap_threshold']:
                    # Gap threshold met, call gap
                    consensus_seq += '-'
                elif not valid_bases:
                    # All positions are N's or gaps but below threshold, call N
                    consensus_seq += 'N'
                else:
                    # Count bases
                    base_counts = Counter(valid_bases)
                    most_common_base, count = base_counts.most_common(1)[0]
                    
                    # Calculate coverage relative to total sequences (including N's from padding)
                    total_sequences = len(bases)
                    coverage = count / total_sequences
                    
                    # Debug: Show what's happening at this position (first 10 positions only)
                    if pos < 10:
                        logger.debug(f"  Pos {pos}: bases={bases}, valid={valid_bases}, counts={dict(base_counts)}, coverage={coverage:.2f}, gap_fraction={gap_fraction:.2f}")
                    
                    # Use more lenient thresholds for Nanopore consensus
                    if coverage >= 0.5:  # 50% of total sequences must agree (more lenient)
                        consensus_seq += most_common_base.upper()
                    elif coverage >= 0.3 and len(base_counts) == 2:
                        # Call ambiguous base for close calls
                        second_most_common_item = base_counts.most_common(2)[1]
                        second_most_common_base = second_most_common_item[0]  # Extract base from (base, count) tuple
                        second_most_common_count = second_most_common_item[1]  # Extract count from (base, count) tuple
                        if second_most_common_count / total_sequences >= 0.2:
                            # Call IUPAC ambiguous base
                            if most_common_base.upper() == 'A' and second_most_common_base.upper() == 'G':
                                consensus_seq += 'R'
                            elif most_common_base.upper() == 'C' and second_most_common_base.upper() == 'T':
                                consensus_seq += 'Y'
                            elif most_common_base.upper() == 'G' and second_most_common_base.upper() == 'C':
                                consensus_seq += 'S'
                            elif most_common_base.upper() == 'A' and second_most_common_base.upper() == 'T':
                                consensus_seq += 'W'
                            elif most_common_base.upper() == 'G' and second_most_common_base.upper() == 'T':
                                consensus_seq += 'K'
                            elif most_common_base.upper() == 'A' and second_most_common_base.upper() == 'C':
                                consensus_seq += 'M'
                            else:
                                consensus_seq += 'N'
                        else:
                            consensus_seq += most_common_base
                    else:
                        # Always call the most common base if available
                        consensus_seq += most_common_base.upper()
            
            # Trim trailing N's
            consensus_seq = consensus_seq.rstrip('N')
            
            if not consensus_seq:
                logger.warning(f"Cluster {cluster_id} consensus is empty after trimming")
                return None
            
            # Create consensus sequence record
            consensus_record = SeqRecord(
                seq=Seq(consensus_seq),
                id=f"cluster_{cluster_id:04d}_nanopore_consensus",
                description=f"Nanopore consensus for cluster {cluster_id} ({len(consensus_seq)} bp, no alignment gaps)"
            )
            
            logger.info(f"Created Nanopore consensus for cluster {cluster_id}: {len(consensus_seq)} bp")
            return consensus_record
            
        except Exception as e:
            logger.error(f"Error creating Nanopore consensus for cluster {cluster_id}: {e}")
            return None
    
    def _calculate_consensus(self, aligned_file, cluster_id):
        """Calculate consensus sequence from aligned reads."""
        try:
            # Read aligned sequences
            alignment = AlignIO.read(aligned_file, 'fasta')
            
            if len(alignment) < 2:
                logger.warning(f"Cluster {cluster_id} has only {len(alignment)} sequence(s), cannot calculate consensus")
                return None
            
            # Calculate consensus using majority voting
            consensus_seq = ""
            alignment_length = alignment.get_alignment_length()
            
            for pos in range(alignment_length):
                # Get all bases at this position
                bases = [str(alignment[i, pos]) for i in range(len(alignment))]
                
                # Count gaps and bases
                gap_count = sum(1 for base in bases if base.upper() == '-')
                gap_fraction = gap_count / len(bases)
                base_counts = Counter([base.upper() for base in bases if base.upper() != '-'])
                
                # Check if gap threshold is met
                if gap_fraction >= self.polishing_params['gap_threshold']:
                    # Gap threshold met, call gap
                    consensus_seq += '-'
                elif not base_counts:
                    # All gaps at this position but below threshold, call N
                    consensus_seq += 'N'
                else:
                    # Find most common base
                    most_common_base, count = base_counts.most_common(1)[0]
                    coverage = count / len(bases)
                    
                    if coverage >= self.polishing_params['min_coverage_threshold']:
                        # Check if we should call an ambiguous base
                        if len(base_counts) > 1:
                            second_most_common = base_counts.most_common(2)[1][1]
                            if second_most_common / len(bases) >= self.polishing_params['ambiguous_threshold']:
                                # Call ambiguous base
                                if most_common_base == 'A' and second_most_common == 'G':
                                    consensus_seq += 'R'
                                elif most_common_base == 'C' and second_most_common == 'T':
                                    consensus_seq += 'Y'
                                elif most_common_base == 'G' and second_most_common == 'C':
                                    consensus_seq += 'S'
                                elif most_common_base == 'A' and second_most_common == 'T':
                                    consensus_seq += 'W'
                                elif most_common_base == 'G' and second_most_common == 'T':
                                    consensus_seq += 'K'
                                elif most_common_base == 'A' and second_most_common == 'C':
                                    consensus_seq += 'M'
                                else:
                                    consensus_seq += 'N'
                            else:
                                consensus_seq += most_common_base
                        else:
                            consensus_seq += most_common_base
                    else:
                        # Coverage too low, call N instead of gap
                        consensus_seq += 'N'
            
            # Remove leading and trailing gaps
            consensus_seq = consensus_seq.strip('-')
            
            if not consensus_seq:
                logger.warning(f"Cluster {cluster_id} consensus is empty after gap trimming")
                return None
            
            # Create consensus sequence record
            consensus_record = SeqRecord(
                seq=Seq(consensus_seq),
                id=f"cluster_{cluster_id:04d}_consensus",
                description=f"Consensus sequence for cluster {cluster_id} ({len(consensus_seq)} bp)"
            )
            
            logger.info(f"Calculated consensus for cluster {cluster_id}: {len(consensus_seq)} bp")
            return consensus_record
            
        except Exception as e:
            logger.error(f"Error calculating consensus for cluster {cluster_id}: {e}")
            return None
    
    def _polish_cluster(self, cluster_file, cluster_id, temp_dir):
        """Polish a cluster by calculating its consensus sequence."""
        if not self.polish_clusters:
            return None
        
        logger.info(f"Polishing cluster {cluster_id} with consensus calculation")
        
        # Debug: Check cluster file content
        try:
            reads = list(SeqIO.parse(cluster_file, 'fasta'))
            logger.info(f"Cluster {cluster_id} contains {len(reads)} sequences")
            for i, read in enumerate(reads):
                logger.debug(f"  Read {i+1}: {read.id} ({len(read.seq)} bp)")
        except Exception as e:
            logger.error(f"Error reading cluster file {cluster_file}: {e}")
            raise RuntimeError(f"Error reading cluster file {cluster_file}: {e}")
        
                # Use MAFFT for alignment (superior to MUSCLE for similar sequences)
        try:
            logger.info(f"Running MAFFT alignment with Nanopore settings for cluster {cluster_id}")
            aligned_file = self._run_mafft_nanopore(cluster_file, temp_dir, cluster_id)
            
            if aligned_file:
                # Calculate consensus from MAFFT alignment
                consensus = self._calculate_consensus(aligned_file, cluster_id)
                if consensus:
                    self.consensus_list.append(consensus)
                    self.cluster_stats['clusters_polished'] += 1
                    
                    # Write consensus to cluster-specific file
                    consensus_file = self.output_dir / f"cluster_{cluster_id:04d}_consensus.fasta"
                    SeqIO.write([consensus], consensus_file, 'fasta')
                    
                    logger.info(f"MAFFT consensus written to: {consensus_file}")
                    return consensus
                else:
                    logger.warning(f"Failed to calculate consensus from MAFFT alignment for cluster {cluster_id}")
            else:
                logger.warning(f"MAFFT alignment failed for cluster {cluster_id}")
        except Exception as e:
            logger.error(f"MAFFT alignment failed for cluster {cluster_id}: {e}")
        
        # Fallback: Create simple consensus without alignment
        logger.info(f"Using fallback simple consensus method for cluster {cluster_id}")
        consensus = self._create_nanopore_consensus(cluster_file, cluster_id)
        if consensus:
            self.consensus_list.append(consensus)
            self.cluster_stats['clusters_polished'] += 1
            
            # Write consensus to cluster-specific file
            consensus_file = self.output_dir / f"cluster_{cluster_id:04d}_consensus.fasta"
            SeqIO.write([consensus], consensus_file, 'fasta')
            
            logger.info(f"Simple consensus written to: {consensus_file}")
            return consensus
        else:
            logger.error(f"Failed to create consensus for cluster {cluster_id} using all methods")
            raise RuntimeError(f"Failed to create consensus for cluster {cluster_id} using all methods")
    
    def cluster_reads(self):
        """Main clustering function using greedy algorithm."""
        logger.info("Starting read clustering process")
        
        # Create temp directory
        temp_dir = self.output_dir / "tmp"
        temp_dir.mkdir(exist_ok=True)
        
        try:
            # Step 1: Sort reads by length (longest first)
            logger.info("Step 1: Sorting reads by length")
            reads = self._sort_reads_by_length(self.input_file)
            
            if not reads:
                logger.warning("No reads found in input file")
                return {}
            
            # Step 2: Choose clustering method
            if self.use_fast_clustering:
                logger.info("Step 2: Using VSEARCH for fast clustering")
                clusters = self._fast_clustering_vsearch(reads, temp_dir)
                cluster_files = self._process_vsearch_clusters(clusters, temp_dir)
            else:
                logger.info("Step 2: Using BLAST for greedy clustering")
                cluster_files = self._blast_based_clustering(reads, temp_dir)
            
            # Step 3: Handle remaining unclustered reads
            if 'unclustered' not in cluster_files:
                unclustered_file = self.output_dir / "unclustered_reads.fasta"
                SeqIO.write(reads, unclustered_file, 'fasta')
                cluster_files["unclustered"] = unclustered_file
                
                self.cluster_stats['reads_unclustered'] = len(reads)
                logger.info(f"Wrote {len(reads)} unclustered reads to: {unclustered_file}")
            
            # Update final statistics
            self.cluster_stats['total_clusters'] = len(cluster_files)
            
            # Write cluster list file
            self._write_cluster_list()
            
            # Write consensus file if polishing was enabled
            if self.polish_clusters:
                self._write_consensus_file()
            
            # Generate summary report
            self._write_summary_report(cluster_files)
            
            logger.info("Clustering completed successfully")
            return cluster_files
            
        except Exception as e:
            logger.error(f"Error during clustering: {e}")
            raise
        finally:
            # Clean up temp directory
            if not self.save_intermediate:
                shutil.rmtree(temp_dir, ignore_errors=True)
                logger.info("Cleaned up temporary files")
    
    def _blast_based_clustering(self, reads, temp_dir):
        """Original BLAST-based greedy clustering algorithm."""
        cluster_id = 1
        cluster_files = {}
        
        while reads and (self.max_clusters is None or cluster_id <= self.max_clusters):
            # Pick the longest remaining read as representative
            representative = reads[0]
            representative_length = len(representative.seq)
            
            logger.info(f"Cluster {cluster_id}: Using {representative.id} ({representative_length} bp) as representative")
            
            # Create temporary representative file
            temp_rep_file = temp_dir / f"representative_{cluster_id}.fasta"
            SeqIO.write([representative], temp_rep_file, 'fasta')
            
            # Create BLAST database for representative
            db_name = f"rep_db_{cluster_id}"
            db_path = self._create_blast_db(temp_rep_file, db_name, temp_dir)
            
            # Create temporary reads file (excluding representative)
            remaining_reads_file = temp_dir / f"remaining_reads_{cluster_id}.fasta"
            SeqIO.write(reads[1:], remaining_reads_file, 'fasta')
            
            # Run BLAST search
            blast_output = temp_dir / f"blast_output_{cluster_id}.txt"
            self._run_blast(remaining_reads_file, db_path, blast_output, cluster_id)
            
            # Parse BLAST results
            hits = self._parse_blast_results(blast_output, representative_length)
            
            # Get reads that meet clustering criteria
            clustered_read_ids = set(hits.keys())
            clustered_reads = [read for read in reads[1:] if read.id in clustered_read_ids]
            
            # Add representative to cluster if it meets minimum size requirement
            total_cluster_size = len(clustered_reads) + 1  # +1 for representative
            
            if total_cluster_size >= self.min_cluster_size:
                # Create cluster file
                cluster_filepath = self._create_cluster(representative, clustered_reads, cluster_id)
                cluster_files[f"cluster_{cluster_id:04d}"] = cluster_filepath
                
                # Polish cluster if requested
                if self.polish_clusters:
                    logger.info(f"Starting polishing for cluster {cluster_id}")
                    polish_result = self._polish_cluster(cluster_filepath, cluster_id, temp_dir)
                    if polish_result:
                        logger.info(f"Successfully polished cluster {cluster_id}")
                    else:
                        logger.warning(f"Failed to polish cluster {cluster_id}")
                else:
                    logger.debug(f"Skipping polishing for cluster {cluster_id} (not enabled)")
                
                # Update statistics
                self.cluster_stats['clusters_created'] += 1
                self.cluster_stats['reads_clustered'] += total_cluster_size
                
                # Remove clustered reads from remaining reads
                reads = self._remove_clustered_reads(reads, clustered_read_ids | {representative.id})
                
                logger.info(f"Cluster {cluster_id} completed: {total_cluster_size} reads")
            else:
                # Skip this representative if cluster is too small
                logger.info(f"Cluster {cluster_id} too small ({total_cluster_size} reads), skipping")
                reads = reads[1:]  # Remove only the representative
            
            cluster_id += 1
        
        return cluster_files
    
    def _process_vsearch_clusters(self, clusters, temp_dir):
        """Process VSEARCH clusters and create cluster files."""
        cluster_files = {}
        cluster_id = 1
        
        for cluster in clusters:
            if len(cluster) >= self.min_cluster_size:
                # Use first read as representative (they're already sorted by length)
                representative = cluster[0]
                cluster_members = cluster[1:]
                
                # Create cluster file
                cluster_filepath = self._create_cluster(representative, cluster_members, cluster_id)
                cluster_files[f"cluster_{cluster_id:04d}"] = cluster_filepath
                
                # Polish cluster if requested
                if self.polish_clusters:
                    logger.info(f"Starting polishing for cluster {cluster_id}")
                    polish_result = self._polish_cluster(cluster_filepath, cluster_id, temp_dir)
                    if polish_result:
                        logger.info(f"Successfully polished cluster {cluster_id}")
                    else:
                        logger.warning(f"Failed to polish cluster {cluster_id}")
                else:
                    logger.debug(f"Skipping polishing for cluster {cluster_id} (not enabled)")
                
                # Update statistics
                self.cluster_stats['clusters_created'] += 1
                self.cluster_stats['reads_clustered'] += len(cluster)
                
                logger.info(f"VSEARCH cluster {cluster_id} completed: {len(cluster)} reads")
                cluster_id += 1
        
        return cluster_files
    
    def _write_consensus_file(self):
        """Write a FASTA file containing all cluster consensus sequences."""
        if not self.polish_clusters or not self.consensus_list:
            return
        
        SeqIO.write(self.consensus_list, self.consensus_file, 'fasta')
        logger.info(f"Consensus sequences written to: {self.consensus_file} ({len(self.consensus_list)} sequences)")
    
    def _write_summary_report(self, cluster_files):
        """Write a summary report of the clustering results."""
        report_file = self.output_dir / "clustering_summary.txt"
        
        with open(report_file, 'w') as f:
            f.write("Read Clustering by BLAST Summary\n")
            f.write("=" * 40 + "\n\n")
            
            # Input information
            f.write("Input Information:\n")
            f.write(f"  Input file: {self.input_file}\n")
            f.write(f"  Output directory: {self.output_dir}\n\n")
            
            # Clustering parameters
            f.write("Clustering Parameters:\n")
            f.write(f"  Clustering method: {'VSEARCH (fast)' if self.use_fast_clustering else 'BLAST (greedy)'}\n")
            f.write(f"  Coverage threshold: {self.coverage_threshold}%\n")
            f.write(f"  Identity threshold: {self.identity_threshold}%\n")
            f.write(f"  Minimum cluster size: {self.min_cluster_size}\n")
            f.write(f"  Maximum clusters: {self.max_clusters if self.max_clusters else 'Unlimited'}\n")
            f.write(f"  Maximum read length: {self.max_read_length if self.max_read_length else 'Unlimited'}\n")
            f.write(f"  Consensus polishing: {'Enabled' if self.polish_clusters else 'Disabled'}\n")
            
            # VSEARCH-specific parameters
            if self.use_fast_clustering:
                f.write("\nVSEARCH Clustering Parameters:\n")
                f.write(f"  Algorithm: {VSEARCH_PARAMS['clustering_algorithm']}\n")
                f.write(f"  Identity threshold: {VSEARCH_PARAMS['identity_threshold']*100:.0f}%\n")
                f.write(f"  Query coverage: {VSEARCH_PARAMS['query_coverage']*100:.0f}%\n")
                f.write(f"  Target coverage: {VSEARCH_PARAMS['target_coverage']*100:.0f}%\n")
                f.write(f"  Strand handling: {VSEARCH_PARAMS['strand']}\n")
                f.write(f"  Size annotations: {'Enabled' if VSEARCH_PARAMS['sizein'] else 'Disabled'}\n")
            
            f.write("\n")
            
            # Polishing parameters (if enabled)
            if self.polish_clusters:
                f.write("Polishing Parameters:\n")
                f.write(f"  Method: MAFFT alignment + Nanopore-optimized consensus\n")
                f.write(f"  Alignment tool: MAFFT (superior to MUSCLE for similar sequences)\n")
                f.write(f"  Coverage threshold: 60% for confident calls\n")
                f.write(f"  Ambiguous base threshold: {self.polishing_params['ambiguous_threshold']*100:.0f}% for IUPAC codes\n")
                f.write(f"  Padding: N's instead of gaps for length differences\n")
                f.write(f"  Max reads for MAFFT: {self.polishing_params['max_reads_for_polishing']} (random sampling for large clusters)\n")
                if self.use_fast_clustering:
                    f.write(f"  VSEARCH optimization: 85% identity clustering + quality filtering for better MAFFT alignment\n")
                f.write("\n")
            
            # Results summary
            f.write("Results Summary:\n")
            f.write(f"  Total reads processed: {self.cluster_stats['total_reads']}\n")
            f.write(f"  Clusters created: {self.cluster_stats['clusters_created']}\n")
            f.write(f"  Reads clustered: {self.cluster_stats['reads_clustered']}\n")
            f.write(f"  Reads unclustered: {self.cluster_stats['reads_unclustered']}\n")
            f.write(f"  Clustering efficiency: {(self.cluster_stats['reads_clustered']/self.cluster_stats['total_reads']*100):.1f}%\n")
            if self.polish_clusters:
                f.write(f"  Clusters polished: {self.cluster_stats['clusters_polished']}\n")
            f.write("\n")
            
            # Cluster details
            f.write("Cluster Details:\n")
            for cluster_name, filepath in cluster_files.items():
                if cluster_name != "unclustered":
                    read_count = len(list(SeqIO.parse(filepath, 'fasta')))
                    f.write(f"  {cluster_name}: {read_count} reads -> {filepath}\n")
                    if self.polish_clusters:
                        consensus_file = self.output_dir / f"{cluster_name}_consensus.fasta"
                        if consensus_file.exists():
                            f.write(f"    Consensus: {consensus_file}\n")
                else:
                    f.write(f"  {cluster_name}: {self.cluster_stats['reads_unclustered']} reads -> {filepath}\n")
            
            # Cluster list information
            f.write(f"\nCluster List:\n")
            f.write(f"  cluster_list.fasta: {len(self.cluster_list)} representatives -> {self.output_dir}/cluster_list.fasta\n")
            f.write(f"  cluster_list_growing.fasta: Real-time growing list -> {self.output_dir}/cluster_list_growing.fasta\n")
            
            # Consensus information (if enabled)
            if self.polish_clusters:
                f.write(f"\nConsensus Sequences:\n")
                f.write(f"  cluster_consensus.fasta: {len(self.consensus_list)} consensus sequences -> {self.output_dir}/cluster_consensus.fasta\n")
        
        logger.info(f"Summary report written to: {report_file}")
    
    def _write_cluster_list(self):
        """Write a FASTA file containing all cluster representatives."""
        cluster_list_file = self.output_dir / "cluster_list.fasta"
        
        if self.cluster_list:
            SeqIO.write(self.cluster_list, cluster_list_file, 'fasta')
            logger.info(f"Cluster list written to: {cluster_list_file} ({len(self.cluster_list)} representatives)")
        else:
            logger.warning("No clusters created, cluster_list.fasta will be empty")
            # Create empty file
            with open(cluster_list_file, 'w') as f:
                pass
        
        # Note: cluster_list_growing.fasta is already complete and up-to-date
        logger.info(f"Growing cluster list available at: {self.growing_cluster_file}")


def main():
    """Main function with command line argument parsing."""
    parser = argparse.ArgumentParser(
        description="Cluster reads using greedy BLAST-based algorithm",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python read_cluster_by_blast.py input.fasta output_dir
  python read_cluster_by_blast.py input.fasta output_dir --coverage 90 --identity 85
  python read_cluster_by_blast.py input.fasta output_dir --min-cluster-size 5 --max-clusters 100
  python read_cluster_by_blast.py input.fasta output_dir --max-read-length 5000
  python read_cluster_by_blast.py input.fasta output_dir --verbose --save-blast-outputs
  python read_cluster_by_blast.py input.fasta output_dir -p --verbose  # MAFFT + Nanopore consensus
  python read_cluster_by_blast.py input.fasta output_dir -f  # Fast clustering with VSEARCH
  python read_cluster_by_blast.py input.fasta output_dir -f -p  # Fast clustering + MAFFT consensus
        """
    )
    
    parser.add_argument('input_file', help='Input FASTA file')
    parser.add_argument('output_dir', help='Output directory for cluster files')
    parser.add_argument('--coverage', type=int, default=95, help='Minimum coverage percentage for clustering (default: 95)')
    parser.add_argument('--identity', type=int, default=90, help='Minimum identity percentage for clustering (default: 90)')
    parser.add_argument('--min-cluster-size', type=int, default=1, help='Minimum number of reads in a cluster (default: 1)')
    parser.add_argument('--max-clusters', type=int, help='Maximum number of clusters to create (default: unlimited)')
    parser.add_argument('--max-read-length', type=int, help='Maximum read length to consider for clustering (default: unlimited)')
    parser.add_argument('-p', '--polish', action='store_true', help='Calculate consensus sequences for clusters using MAFFT alignment + Nanopore-optimized method')
    parser.add_argument('-f', '--fast', action='store_true', help='Use VSEARCH for fast clustering instead of BLAST (10-100x faster, optimized for MAFFT consensus)')
    parser.add_argument('--debug', action='store_true', help='Enable debug logging')
    parser.add_argument('--verbose', action='store_true', help='Print detailed output')
    parser.add_argument('--save-blast-outputs', action='store_true', help='Save BLAST output files for debugging')
    parser.add_argument('--save-intermediate', action='store_true', help='Save intermediate files')
    
    args = parser.parse_args()
    
    # Set up logging level
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
        logger.setLevel(logging.DEBUG)
    
    try:
        # Create clusterer instance
        clusterer = ReadClusterer(
            input_file=args.input_file,
            output_dir=args.output_dir,
            coverage_threshold=args.coverage,
            identity_threshold=args.identity,
            min_cluster_size=args.min_cluster_size,
            max_clusters=args.max_clusters,
            max_read_length=args.max_read_length,
            save_blast_outputs=args.save_blast_outputs,
            verbose=args.verbose,
            save_intermediate=args.save_intermediate,
            polish_clusters=args.polish,
            use_fast_clustering=args.fast
        )
        
        # Run clustering
        cluster_files = clusterer.cluster_reads()
        
        logger.info("Clustering completed successfully!")
        logger.info(f"Cluster files created in: {args.output_dir}")
        
        if args.polish:
            logger.info("Consensus sequences calculated using MAFFT and saved to cluster_consensus.fasta")
        
    except Exception as e:
        logger.error(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
