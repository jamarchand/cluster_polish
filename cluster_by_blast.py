#!/usr/bin/env python3
"""
BLAST-based clustering module.

This module implements the greedy BLAST-based clustering algorithm that groups reads
based on high sequence similarity using BLAST. It picks the longest read as a 
representative, finds all reads that cover >95% of it, and creates clusters.
"""

import os
import subprocess
import shutil
import logging
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from params import BLAST_PARAMS

logger = logging.getLogger(__name__)


class BlastClusterer:
    """BLAST-based clustering implementation."""
    
    def __init__(self, blast_params=None):
        """
        Initialize BLAST clusterer.
        
        Args:
            blast_params (dict): BLAST parameters (uses defaults if None)
        """
        self.blast_params = blast_params or BLAST_PARAMS.copy()
    
    def create_blast_db(self, query_file, db_name, temp_dir):
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
    
    def run_blast(self, query_file, db_path, output_file, cluster_id=None, save_blast_outputs=False, output_dir=None):
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
            if save_blast_outputs and cluster_id and output_dir:
                blast_outputs_dir = output_dir / "blast_outputs"
                blast_outputs_dir.mkdir(exist_ok=True)
                
                local_blast_output = blast_outputs_dir / f"blast_cluster_{cluster_id}_output.txt"
                shutil.copy2(output_file, local_blast_output)
                logger.info(f"Saved BLAST output to: {local_blast_output}")
                
        except subprocess.CalledProcessError as e:
            logger.error(f"BLAST search failed: {e}")
            raise
    
    def parse_blast_results(self, blast_output, query_length, identity_threshold, coverage_threshold, verbose=False):
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
                    if pident >= identity_threshold and query_coverage_percent >= coverage_threshold:
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
                        
                        if verbose:
                            logger.debug(f"Hit for read {qseqid}: identity={pident}%, coverage={query_coverage_percent:.1f}%")
        
        logger.info(f"Found {len(hits)} reads meeting clustering criteria")
        return hits
    
    def cluster_reads(self, reads, temp_dir, output_dir, coverage_threshold, identity_threshold, 
                     min_cluster_size, max_clusters, save_blast_outputs=False, verbose=False, output_longest_reads=False):
        """
        Perform BLAST-based greedy clustering.
        
        Args:
            reads: List of SeqRecord objects sorted by length (longest first)
            temp_dir: Temporary directory for intermediate files
            output_dir: Output directory for cluster files
            coverage_threshold: Minimum coverage percentage for clustering
            identity_threshold: Minimum identity percentage for clustering
            min_cluster_size: Minimum number of reads in a cluster
            max_clusters: Maximum number of clusters to create
            save_blast_outputs: Whether to save BLAST output files
            verbose: Whether to print detailed output
            output_longest_reads: Whether to collect longest reads from clusters
            
        Returns:
            tuple: (cluster_files_dict, cluster_list, longest_reads_list)
        """
        cluster_id = 1
        cluster_files = {}
        cluster_list = []
        longest_reads_list = []
        growing_cluster_file = output_dir / "cluster_list_growing.fasta"
        
        # Initialize growing cluster list file
        if growing_cluster_file.exists():
            growing_cluster_file.unlink()
        
        while reads and (max_clusters is None or cluster_id <= max_clusters):
            # Pick the longest remaining read as representative
            representative = reads[0]
            representative_length = len(representative.seq)
            
            logger.info(f"Cluster {cluster_id}: Using {representative.id} ({representative_length} bp) as representative")
            
            # Create temporary representative file
            temp_rep_file = temp_dir / f"representative_{cluster_id}.fasta"
            SeqIO.write([representative], temp_rep_file, 'fasta')
            
            # Create BLAST database for representative
            db_name = f"rep_db_{cluster_id}"
            db_path = self.create_blast_db(temp_rep_file, db_name, temp_dir)
            
            # Create temporary reads file (excluding representative)
            remaining_reads_file = temp_dir / f"remaining_reads_{cluster_id}.fasta"
            SeqIO.write(reads[1:], remaining_reads_file, 'fasta')
            
            # Run BLAST search
            blast_output = temp_dir / f"blast_output_{cluster_id}.txt"
            self.run_blast(remaining_reads_file, db_path, blast_output, cluster_id, save_blast_outputs, output_dir)
            
            # Parse BLAST results
            hits = self.parse_blast_results(blast_output, representative_length, identity_threshold, coverage_threshold, verbose)
            
            # Get reads that meet clustering criteria
            clustered_read_ids = set(hits.keys())
            clustered_reads = [read for read in reads[1:] if read.id in clustered_read_ids]
            
            # Add representative to cluster if it meets minimum size requirement
            total_cluster_size = len(clustered_reads) + 1  # +1 for representative
            
            if total_cluster_size >= min_cluster_size:
                # Create cluster file
                cluster_filepath = self._create_cluster(representative, clustered_reads, cluster_id, output_dir, cluster_list, growing_cluster_file)
                cluster_files[f"cluster_{cluster_id:04d}"] = cluster_filepath
                
                # Collect longest read if requested
                if output_longest_reads:
                    longest_read = max([representative] + clustered_reads, key=lambda x: len(x.seq))
                    longest_read_copy = SeqRecord(
                        seq=longest_read.seq,
                        id=f"cluster_{cluster_id:04d}_longest",
                        description=f"Longest read from cluster {cluster_id} (original: {longest_read.id}, {len(longest_read.seq)} bp)"
                    )
                    longest_reads_list.append(longest_read_copy)
                
                logger.info(f"Cluster {cluster_id} completed: {total_cluster_size} reads")
                
                # Remove clustered reads from remaining reads
                reads = self._remove_clustered_reads(reads, clustered_read_ids | {representative.id})
            else:
                # Skip this representative if cluster is too small
                logger.info(f"Cluster {cluster_id} too small ({total_cluster_size} reads), skipping")
                reads = reads[1:]  # Remove only the representative
            
            cluster_id += 1
        
        return cluster_files, cluster_list, longest_reads_list
    
    def _create_cluster(self, representative_read, cluster_reads, cluster_id, output_dir, cluster_list, growing_cluster_file):
        """Create a cluster file with the representative and cluster reads."""
        cluster_filename = f"cluster_{cluster_id:04d}.fasta"
        cluster_filepath = output_dir / cluster_filename
        
        # Write representative read first, then cluster reads
        all_reads = [representative_read] + cluster_reads
        SeqIO.write(all_reads, cluster_filepath, 'fasta')
        
        # Add representative to cluster list for tracking
        representative_copy = SeqRecord(
            seq=representative_read.seq,
            id=f"cluster_{cluster_id:04d}_representative",
            description=f"Representative for cluster {cluster_id} (original: {representative_read.id})"
        )
        cluster_list.append(representative_copy)
        
        # Write to growing cluster file immediately
        with open(growing_cluster_file, 'a') as f:
            SeqIO.write(representative_copy, f, 'fasta')
        
        logger.info(f"Added cluster {cluster_id} representative to growing list: {growing_cluster_file}")
        logger.info(f"Created cluster {cluster_id}: {len(all_reads)} reads -> {cluster_filepath}")
        
        return cluster_filepath
    
    def _remove_clustered_reads(self, reads, clustered_read_ids):
        """Remove clustered reads from the reads list."""
        remaining_reads = [read for read in reads if read.id not in clustered_read_ids]
        removed_count = len(reads) - len(remaining_reads)
        
        logger.info(f"Removed {removed_count} clustered reads, {len(remaining_reads)} remaining")
        
        return remaining_reads
