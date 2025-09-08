#!/usr/bin/env python3
"""
MAFFT-based polishing module.

This module implements consensus sequence calculation using MAFFT alignment
and Nanopore-optimized consensus methods for error correction.
"""

import os
import subprocess
import shutil
import random
import logging
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import AlignIO
from collections import Counter

from params import POLISHING_PARAMS

logger = logging.getLogger(__name__)


class MafftPolisher:
    """MAFFT-based consensus polishing implementation."""
    
    def __init__(self, polishing_params=None):
        """
        Initialize MAFFT polisher.
        
        Args:
            polishing_params (dict): Polishing parameters (uses defaults if None)
        """
        self.polishing_params = polishing_params or POLISHING_PARAMS.copy()
    
    def check_mafft_availability(self):
        """Check if MAFFT is available and get version."""
        try:
            result = subprocess.run(['mafft', '--version'], capture_output=True, text=True, check=True)
            logger.info(f"MAFFT version: {result.stdout.strip()}")
            return True
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            logger.warning(f"MAFFT not available: {e}")
            return False
    
    def run_mafft_nanopore(self, cluster_file, temp_dir, cluster_id):
        """Run MAFFT with Nanopore-optimized settings for error correction."""
        if not self.check_mafft_availability():
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
    
    def calculate_consensus(self, aligned_file, cluster_id):
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
    
    def create_nanopore_consensus(self, cluster_file, cluster_id):
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
    
    def polish_cluster(self, cluster_file, cluster_id, temp_dir, output_dir):
        """
        Polish a cluster by calculating its consensus sequence.
        
        Args:
            cluster_file: Path to cluster FASTA file
            cluster_id: Cluster identifier
            temp_dir: Temporary directory for intermediate files
            output_dir: Output directory for consensus files
            
        Returns:
            SeqRecord: Consensus sequence record, or None if failed
        """
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
        
        # Try MAFFT alignment first
        try:
            logger.info(f"Running MAFFT alignment with Nanopore settings for cluster {cluster_id}")
            aligned_file = self.run_mafft_nanopore(cluster_file, temp_dir, cluster_id)
            
            if aligned_file:
                # Calculate consensus from MAFFT alignment
                consensus = self.calculate_consensus(aligned_file, cluster_id)
                if consensus:
                    # Write consensus to cluster-specific file
                    consensus_file = output_dir / f"cluster_{cluster_id:04d}_consensus.fasta"
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
        consensus = self.create_nanopore_consensus(cluster_file, cluster_id)
        if consensus:
            # Write consensus to cluster-specific file
            consensus_file = output_dir / f"cluster_{cluster_id:04d}_consensus.fasta"
            SeqIO.write([consensus], consensus_file, 'fasta')
            
            logger.info(f"Simple consensus written to: {consensus_file}")
            return consensus
        else:
            logger.error(f"Failed to create consensus for cluster {cluster_id} using all methods")
            raise RuntimeError(f"Failed to create consensus for cluster {cluster_id} using all methods")
