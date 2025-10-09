#!/usr/bin/env python3
"""
Cluster Polish Tool

This script implements a modular clustering pipeline that groups reads based on 
high sequence similarity using VSEARCH, with optional RACON-based consensus polishing.

Author: Assistant
Date: 2024
"""

import os
import sys
import shutil
import subprocess
import tempfile
from pathlib import Path
from Bio import SeqIO
import argparse
import logging

# Import modular components
from params import (
    VSEARCH_PARAMS, RACON_PARAMS, DEBUG_PARAMS,
    REDUNDANT_REMOVAL_PARAMS, DEFAULT_CLUSTERING_ALGORITHM, DEFAULT_POLISHING_ALGORITHM
)
from cluster_by_vsearch import VsearchClusterer
from polish_by_racon import RaconPolisher

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class ReadClusterer:
    def __init__(self, input_file, output_dir, coverage_threshold=None, identity_threshold=None,
                 min_cluster_size=None, max_clusters=None, max_read_length=None,
                 verbose=None, save_intermediate=None, polish_clusters=None,
                 clustering_algorithm=None, polishing_algorithm=None, output_longest_reads=None,
                 remove_redundant_consensus=None):
        """
        Initialize the read clusterer with input file and parameters.
        
        Args:
            input_file (str): Path to input FASTA file
            output_dir (str): Directory for output files
            coverage_threshold (int): Minimum coverage percentage for clustering
            identity_threshold (int): Minimum identity percentage for clustering
            min_cluster_size (int): Minimum reads in a cluster
            max_clusters (int): Maximum number of clusters
            verbose (bool): Print detailed output
            save_intermediate (bool): Save intermediate files
            polish_clusters (bool): Whether to calculate consensus sequences for clusters
            clustering_algorithm (str): Clustering algorithm to use ('vsearch')
            polishing_algorithm (str): Polishing algorithm to use ('racon')
            output_longest_reads (bool): Whether to output longest read from each cluster
            remove_redundant_consensus (bool): Whether to remove redundant consensus sequences
        """
        self.input_file = Path(input_file)
        self.output_dir = Path(output_dir)
        
        # Use provided parameters or defaults
        # Fallbacks mirror argparse defaults
        self.coverage_threshold = coverage_threshold if coverage_threshold is not None else 95
        self.identity_threshold = identity_threshold if identity_threshold is not None else 90
        self.min_cluster_size = min_cluster_size if min_cluster_size is not None else 1
        self.max_clusters = max_clusters if max_clusters is not None else None
        self.max_read_length = max_read_length if max_read_length is not None else None
        self.verbose = verbose if verbose is not None else DEBUG_PARAMS['verbose']
        self.save_intermediate = save_intermediate if save_intermediate is not None else DEBUG_PARAMS['save_intermediate']
        self.polish_clusters = polish_clusters if polish_clusters is not None else False
        self.output_longest_reads = output_longest_reads if output_longest_reads is not None else DEBUG_PARAMS['output_longest_reads']
        self.remove_redundant_consensus = remove_redundant_consensus if remove_redundant_consensus is not None else REDUNDANT_REMOVAL_PARAMS['enabled']
        
        # Algorithm selection
        self.clustering_algorithm = clustering_algorithm or DEFAULT_CLUSTERING_ALGORITHM
        self.polishing_algorithm = polishing_algorithm or DEFAULT_POLISHING_ALGORITHM
        
        # Initialize algorithm-specific components
        self._initialize_components()
        
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
            'clusters_polished': 0,
            'consensus_sequences': 0,
            'redundant_removed': 0,
            'final_consensus_sequences': 0
        }
        
        # Initialize cluster list for tracking parent sequences
        self.cluster_list = []
        
        # Initialize growing cluster list file
        self.growing_cluster_file = self.output_dir / "cluster_list_growing.fasta"
        
        # Initialize longest reads list if requested
        if self.output_longest_reads:
            self.longest_reads_list = []
            self.longest_reads_file = self.output_dir / "cluster_longest_reads.fasta"
        
        # Initialize polished consensus file
        if self.polish_clusters:
            self.consensus_file = self.output_dir / "cluster_consensus.fasta"
            self.consensus_list = []
    
    def _initialize_components(self):
        """Initialize clustering and polishing components based on selected algorithms."""
        # Initialize clustering component (VSEARCH only)
        if self.clustering_algorithm == 'vsearch':
            self.clusterer = VsearchClusterer(VSEARCH_PARAMS)
        else:
            raise ValueError(f"Unknown clustering algorithm: {self.clustering_algorithm}")
        
        # Initialize polishing component
        if self.polish_clusters:
            if self.polishing_algorithm == 'racon':
                self.polisher = RaconPolisher(RACON_PARAMS)
            else:
                raise ValueError(f"Unknown polishing algorithm: {self.polishing_algorithm}")
        else:
            self.polisher = None
        
    def _validate_inputs(self):
        """Validate that input file exists and is readable."""
        if not self.input_file.exists():
            raise FileNotFoundError(f"Input file not found: {self.input_file}")
        
        # Check if input is FASTA
        self.input_format = self._detect_format(self.input_file)
        if self.input_format != 'fasta':
            raise ValueError(f"Input file must be in FASTA format, detected: {self.input_format}")
        
        # Note: RACON is required for consensus calculation
        if self.polish_clusters:
            if self.polishing_algorithm == 'racon':
                logger.info("Nanopore consensus calculation enabled (RACON and minimap2 required for alignment)")
        
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
    
    def cluster_reads(self):
        """Main clustering function using modular components."""
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
            
            # Step 2: Run clustering using selected algorithm
            logger.info(f"Step 2: Using VSEARCH for clustering")
            cluster_files, cluster_list, longest_reads_list = self.clusterer.cluster_reads(
                reads, temp_dir, self.output_dir,
                self.min_cluster_size, use_quality_filtering=True, output_longest_reads=self.output_longest_reads
            )
            
            # Update cluster list and longest reads list
            self.cluster_list = cluster_list
            if self.output_longest_reads and hasattr(self, 'longest_reads_list'):
                self.longest_reads_list = longest_reads_list
            
            # Step 3: Polish clusters if requested
            if self.polish_clusters and self.polisher:
                logger.info("Step 3: Polishing clusters with consensus calculation")
                self._polish_all_clusters(cluster_files, temp_dir)
            
            # Step 4: Handle remaining unclustered reads
            if 'unclustered' not in cluster_files:
                unclustered_file = self.output_dir / "unclustered_reads.fasta"
                SeqIO.write(reads, unclustered_file, 'fasta')
                cluster_files["unclustered"] = unclustered_file
                
                self.cluster_stats['reads_unclustered'] = len(reads)
                logger.info(f"Wrote {len(reads)} unclustered reads to: {unclustered_file}")
            
            # Update final statistics
            self.cluster_stats['total_clusters'] = len(cluster_files)
            self.cluster_stats['clusters_created'] = len([k for k in cluster_files.keys() if k.startswith('cluster_')])
            self.cluster_stats['reads_clustered'] = sum(
                len(list(SeqIO.parse(f, 'fasta'))) 
                for k, f in cluster_files.items() 
                if k.startswith('cluster_')
            )
            
            # Write cluster list file
            self._write_cluster_list()
            
            # Write consensus file if polishing was enabled
            if self.polish_clusters:
                self._write_consensus_file()
                
                # Remove redundant consensus sequences if enabled
                if self.remove_redundant_consensus:
                    self._remove_redundant_consensus()
            
            # Write longest reads file if requested
            if self.output_longest_reads:
                self._write_longest_reads_file()
            
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
    
    def _polish_all_clusters(self, cluster_files, temp_dir):
        """Polish all clusters using the selected polishing algorithm."""
        for cluster_name, cluster_file in cluster_files.items():
            if cluster_name.startswith('cluster_'):
                # Extract cluster ID from cluster name
                cluster_id = int(cluster_name.split('_')[1])
                
                logger.info(f"Starting polishing for cluster {cluster_id}")
                try:
                    polish_result = self.polisher.polish_cluster(cluster_file, cluster_id, temp_dir, self.output_dir)
                    if polish_result:
                        self.consensus_list.append(polish_result)
                        self.cluster_stats['clusters_polished'] += 1
                        logger.info(f"Successfully polished cluster {cluster_id}")
                    else:
                        logger.warning(f"Failed to polish cluster {cluster_id}")
                except Exception as e:
                    logger.error(f"Error polishing cluster {cluster_id}: {e}")
    
    def _write_consensus_file(self):
        """Write a FASTA file containing all cluster consensus sequences."""
        if not self.polish_clusters or not self.consensus_list:
            return
        
        SeqIO.write(self.consensus_list, self.consensus_file, 'fasta')
        self.cluster_stats['consensus_sequences'] = len(self.consensus_list)
        logger.info(f"Consensus sequences written to: {self.consensus_file} ({len(self.consensus_list)} sequences)")
    
    def _remove_redundant_consensus(self):
        """Remove redundant consensus sequences using VSEARCH."""
        if not self.remove_redundant_consensus or not self.consensus_list:
            return
        
        logger.info("Starting redundant consensus sequence removal using VSEARCH")
        
        # Create temporary files for VSEARCH
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_dir = Path(temp_dir)
            input_file = temp_dir / "consensus_input.fasta"
            output_file = temp_dir / "consensus_deduplicated.fasta"
            uc_file = temp_dir / "consensus_clusters.uc"
            
            # Write consensus sequences to temporary file
            SeqIO.write(self.consensus_list, input_file, 'fasta')
            
            # Build VSEARCH command for clustering
            vsearch_cmd = [
                'vsearch',
                '--cluster_fast', str(input_file),
                '--id', str(REDUNDANT_REMOVAL_PARAMS['identity_threshold']),
                '--query_cov', str(REDUNDANT_REMOVAL_PARAMS['query_coverage']),
                '--target_cov', str(REDUNDANT_REMOVAL_PARAMS['target_coverage']),
                '--strand', REDUNDANT_REMOVAL_PARAMS['strand'],
                '--centroids', str(output_file),
                '--uc', str(uc_file),
                '--threads', '1'
            ]
            
            try:
                # Run VSEARCH
                result = subprocess.run(vsearch_cmd, capture_output=True, text=True, check=True)
                
                if self.verbose:
                    logger.debug(f"VSEARCH output: {result.stdout}")
                    if result.stderr:
                        logger.debug(f"VSEARCH stderr: {result.stderr}")
                
                # Read deduplicated sequences
                deduplicated_sequences = list(SeqIO.parse(output_file, 'fasta'))
                
                # Update statistics
                original_count = len(self.consensus_list)
                final_count = len(deduplicated_sequences)
                removed_count = original_count - final_count
                
                self.cluster_stats['redundant_removed'] = removed_count
                self.cluster_stats['final_consensus_sequences'] = final_count
                
                # Write deduplicated consensus file
                deduplicated_file = self.output_dir / f"cluster_consensus{REDUNDANT_REMOVAL_PARAMS['output_file_suffix']}.fasta"
                SeqIO.write(deduplicated_sequences, deduplicated_file, 'fasta')
                
                logger.info(f"Redundant consensus removal completed:")
                logger.info(f"  Original consensus sequences: {original_count}")
                logger.info(f"  Redundant sequences removed: {removed_count}")
                logger.info(f"  Final consensus sequences: {final_count}")
                logger.info(f"  Deduplicated file: {deduplicated_file}")
                
                # Update consensus list with deduplicated sequences
                self.consensus_list = deduplicated_sequences
                
            except subprocess.CalledProcessError as e:
                logger.error(f"VSEARCH failed with return code {e.returncode}")
                logger.error(f"VSEARCH stderr: {e.stderr}")
                logger.warning("Continuing without redundant sequence removal")
            except FileNotFoundError:
                logger.error("VSEARCH not found in PATH. Please install VSEARCH to use redundant sequence removal.")
                logger.warning("Continuing without redundant sequence removal")
    
    def _write_longest_reads_file(self):
        """Write a FASTA file containing the longest read from each cluster."""
        if not self.output_longest_reads or not hasattr(self, 'longest_reads_list') or not self.longest_reads_list:
            return
        
        SeqIO.write(self.longest_reads_list, self.longest_reads_file, 'fasta')
        logger.info(f"Longest reads written to: {self.longest_reads_file} ({len(self.longest_reads_list)} sequences)")
    
    def _write_summary_report(self, cluster_files):
        """Write a summary report of the clustering results."""
        report_file = self.output_dir / "clustering_summary.txt"
        
        with open(report_file, 'w') as f:
            f.write("Cluster Polish Summary\n")
            f.write("=" * 40 + "\n\n")
            
            # Input information
            f.write("Input Information:\n")
            f.write(f"  Input file: {self.input_file}\n")
            f.write(f"  Output directory: {self.output_dir}\n\n")
            
            # Clustering parameters
            f.write("Clustering Parameters:\n")
            f.write(f"  Clustering algorithm: {self.clustering_algorithm.upper()}\n")
            f.write(f"  Coverage threshold: {self.coverage_threshold}%\n")
            f.write(f"  Identity threshold: {self.identity_threshold}%\n")
            f.write(f"  Minimum cluster size: {self.min_cluster_size}\n")
            f.write(f"  Maximum clusters: {self.max_clusters if self.max_clusters else 'Unlimited'}\n")
            f.write(f"  Maximum read length: {self.max_read_length if self.max_read_length else 'Unlimited'}\n")
            f.write(f"  Consensus polishing: {'Enabled' if self.polish_clusters else 'Disabled'}\n")
            if self.polish_clusters:
                f.write(f"  Polishing algorithm: {self.polishing_algorithm.upper()}\n")
            f.write(f"  Redundant consensus removal: {'Enabled' if self.remove_redundant_consensus else 'Disabled'}\n")
            
            # VSEARCH-specific parameters
            if self.clustering_algorithm == 'vsearch':
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
                f.write(f"  Method: RACON alignment + Nanopore-optimized consensus\n")
                f.write(f"  Alignment tool: minimap2 + RACON (optimized for long-read consensus)\n")
                f.write(f"  Match score: {RACON_PARAMS['match']}\n")
                f.write(f"  Mismatch penalty: {RACON_PARAMS['mismatch']}\n")
                f.write(f"  Gap penalty: {RACON_PARAMS['gap']}\n")
                f.write(f"  Window length: {RACON_PARAMS['window_length']}\n")
                f.write(f"  Quality threshold: {RACON_PARAMS['quality_threshold']}\n")
                f.write(f"  Error threshold: {RACON_PARAMS['error_threshold']}\n")
                f.write(f"  Max reads for RACON: {RACON_PARAMS['max_reads_for_polishing']} (random sampling for large clusters)\n")
                f.write(f"  VSEARCH optimization: Quality filtering for better alignment\n")
                f.write("\n")
            
            # Redundant removal parameters (if enabled)
            if self.remove_redundant_consensus:
                f.write("Redundant Consensus Removal Parameters:\n")
                f.write(f"  Method: VSEARCH cluster_fast algorithm\n")
                f.write(f"  Identity threshold: {REDUNDANT_REMOVAL_PARAMS['identity_threshold']*100:.0f}%\n")
                f.write(f"  Query coverage: {REDUNDANT_REMOVAL_PARAMS['query_coverage']*100:.0f}%\n")
                f.write(f"  Target coverage: {REDUNDANT_REMOVAL_PARAMS['target_coverage']*100:.0f}%\n")
                f.write(f"  Strand handling: {REDUNDANT_REMOVAL_PARAMS['strand']}\n")
                f.write(f"  Output suffix: {REDUNDANT_REMOVAL_PARAMS['output_file_suffix']}\n")
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
                f.write(f"  Consensus sequences created: {self.cluster_stats['consensus_sequences']}\n")
                if self.remove_redundant_consensus:
                    f.write(f"  Redundant sequences removed: {self.cluster_stats['redundant_removed']}\n")
                    f.write(f"  Final consensus sequences: {self.cluster_stats['final_consensus_sequences']}\n")
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
                f.write(f"  cluster_consensus.fasta: {self.cluster_stats['consensus_sequences']} consensus sequences -> {self.output_dir}/cluster_consensus.fasta\n")
                if self.remove_redundant_consensus:
                    deduplicated_file = self.output_dir / f"cluster_consensus{REDUNDANT_REMOVAL_PARAMS['output_file_suffix']}.fasta"
                    f.write(f"  cluster_consensus_deduplicated.fasta: {self.cluster_stats['final_consensus_sequences']} deduplicated sequences -> {deduplicated_file}\n")
            
            # Longest reads information (if enabled)
            if self.output_longest_reads and hasattr(self, 'longest_reads_list') and self.longest_reads_list:
                f.write(f"\nLongest Reads:\n")
                f.write(f"  cluster_longest_reads.fasta: {len(self.longest_reads_list)} longest reads -> {self.output_dir}/cluster_longest_reads.fasta\n")
        
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


def cluster_polish():
    """Main function with command line argument parsing."""
    parser = argparse.ArgumentParser(
        description="Modular read clustering pipeline with VSEARCH clustering and optional RACON polishing",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python cluster_polish.py input.fasta output_dir
  python cluster_polish.py input.fasta output_dir --coverage 90 --identity 85
  python cluster_polish.py input.fasta output_dir --min-cluster-size 5 --max-clusters 100
  python cluster_polish.py input.fasta output_dir --max-read-length 5000
  python cluster_polish.py input.fasta output_dir -p --verbose  # RACON consensus
  python cluster_polish.py input.fasta output_dir --output-longest-reads  # Output longest read from each cluster
  python cluster_polish.py input.fasta output_dir -p --remove-redundant-consensus  # Polish + remove redundant consensus
        """
    )
    
    parser.add_argument('input_file', help='Input FASTA file')
    parser.add_argument('output_dir', help='Output directory for cluster files')
    parser.add_argument('--coverage', type=int, default=95, help='Minimum coverage percentage for clustering (default: 95)')
    parser.add_argument('--identity', type=int, default=90, help='Minimum identity percentage for clustering (default: 90)')
    parser.add_argument('--min-cluster-size', type=int, default=1, help='Minimum number of reads in a cluster (default: 1)')
    parser.add_argument('--max-clusters', type=int, help='Maximum number of clusters to create (default: unlimited)')
    parser.add_argument('--max-read-length', type=int, help='Maximum read length to consider for clustering (default: unlimited)')
    parser.add_argument('-p', '--polish', action='store_true', help='Calculate consensus sequences for clusters using RACON alignment + Nanopore-optimized method')
    parser.add_argument('--debug', action='store_true', help='Enable debug logging')
    parser.add_argument('--verbose', action='store_true', help='Print detailed output')
    parser.add_argument('--save-intermediate', action='store_true', help='Save intermediate files')
    parser.add_argument('--output-longest-reads', action='store_true', help='Output longest read from each cluster to separate FASTA file')
    parser.add_argument('--remove-redundant-consensus', action='store_true', help='Remove redundant consensus sequences using VSEARCH after polishing')
    
    args = parser.parse_args()
    
    # Set up logging level
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
        logger.setLevel(logging.DEBUG)
    
    try:
        # Clustering algorithm is fixed to VSEARCH
        clustering_algorithm = 'vsearch'
        
        # Create clusterer instance
        clusterer = ReadClusterer(
            input_file=args.input_file,
            output_dir=args.output_dir,
            coverage_threshold=args.coverage,
            identity_threshold=args.identity,
            min_cluster_size=args.min_cluster_size,
            max_clusters=args.max_clusters,
            max_read_length=args.max_read_length,
            verbose=args.verbose,
            save_intermediate=args.save_intermediate,
            polish_clusters=args.polish,
            clustering_algorithm=clustering_algorithm,
            polishing_algorithm='racon',
            output_longest_reads=args.output_longest_reads,
            remove_redundant_consensus=args.remove_redundant_consensus
        )
        
        # Run clustering
        cluster_files = clusterer.cluster_reads()
        
        logger.info("Clustering completed successfully!")
        logger.info(f"Cluster files created in: {args.output_dir}")
        
        if args.polish:
            logger.info(f"Consensus sequences calculated using RACON and saved to cluster_consensus.fasta")
            if args.remove_redundant_consensus:
                logger.info(f"Redundant consensus sequences removed and saved to cluster_consensus_deduplicated.fasta")
        
    except Exception as e:
        logger.error(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    cluster_polish()