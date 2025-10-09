# Cluster Polish Tool

A modular read clustering pipeline that groups reads based on high sequence similarity using VSEARCH, with optional RACON-based consensus polishing and redundant sequence removal.

## Features

- **Clustering**: VSEARCH-based clustering
- **Consensus Polishing**: Optional RACON-based consensus calculation
- **Redundant Removal**: VSEARCH-based deduplication of consensus sequences
- **Flexible Parameters**: Configurable thresholds and filtering options
- **Comprehensive Reporting**: Detailed summary reports and statistics
- **Quality Control**: Longest read extraction and intermediate file options

## Installation

### Prerequisites

- Python 3.7+
- BioPython
- VSEARCH (for clustering and redundant removal)
- RACON + minimap2 (for polishing)

### Install Dependencies

```bash
# Install Python dependencies
pip install biopython

# Install external tools (examples for different systems)
# Ubuntu/Debian
sudo apt-get install vsearch racon minimap2

# macOS with Homebrew
brew install vsearch racon minimap2

# Conda
conda install -c bioconda vsearch racon minimap2
```

## Quick Start

```bash
# Basic clustering (VSEARCH)
python cluster_polish.py input.fasta output_dir

# Clustering with RACON consensus polishing
python cluster_polish.py input.fasta output_dir -p

# Advanced example with more options
python cluster_polish.py input.fasta output_dir --output-longest-reads --min-cluster-size 15 --max-read-length 10000 -p
```

## Usage

### Basic Syntax

```bash
python cluster_polish.py <input_file> <output_dir> [options]
```

### Required Arguments

- `input_file`: Input FASTA file containing reads to cluster
- `output_dir`: Output directory for cluster files

### Optional Arguments

#### Clustering Parameters
- `--coverage`: Minimum coverage percentage for clustering (default: 95)
- `--identity`: Minimum identity percentage for clustering (default: 90)
- `--min-cluster-size`: Minimum number of reads in a cluster (default: 1)
- `--max-clusters`: Maximum number of clusters to create (default: unlimited)
- `--max-read-length`: Maximum read length to consider (default: unlimited)

#### Polishing
- `-p, --polish`: Calculate consensus sequences for clusters (RACON)

#### Advanced Features
- `--remove-redundant-consensus`: Remove redundant consensus sequences using VSEARCH
- `--output-longest-reads`: Output longest read from each cluster
- `--save-intermediate`: Save intermediate files for debugging

#### Output Control
- `--verbose`: Print detailed output
- `--debug`: Enable debug logging

## Examples

### Example 1: Basic Clustering (VSEARCH)
```bash
python cluster_polish.py reads.fasta clusters_basic
```

### Example 2: VSEARCH Clustering with RACON Consensus
```bash
python cluster_polish.py reads.fasta clusters_vsearch -p
```

### Example 3: Advanced Pipeline with Options
```bash
python cluster_polish.py inputs.fasta output_dir_name --output-longest-reads --min-cluster-size 15 --max-read-length 10000 -p
```

### Example 4: High-Quality Clustering with Redundant Removal
```bash
python cluster_polish.py reads.fasta clusters_high_quality --coverage 98 --identity 95 --min-cluster-size 10 -p --remove-redundant-consensus
```

### Example 4: RACON Polishing with Debugging
```bash
python cluster_polish.py reads.fasta clusters_racon --verbose --save-intermediate -p
```

### Example 5: Large Dataset Processing
```bash
python cluster_polish.py large_dataset.fasta clusters_large --max-clusters 1000 --max-read-length 5000 -p --output-longest-reads
```

### Example 7: Conservative Clustering
```bash
python cluster_polish.py inputs/kpc_hq.fasta outputs/clusters_conservative --coverage 99 --identity 98 --min-cluster-size 20 --max-clusters 100 --max-read-length 1000
```

### Example 6: Complete Pipeline with Deduplication
```bash
python cluster_polish.py reads.fasta final_clusters --remove-redundant-consensus --output-longest-reads --verbose -p
```

## Output Files

### Core Output Files
- `cluster_*.fasta`: Individual cluster files
- `cluster_list.fasta`: Representative sequences from each cluster
- `cluster_list_growing.fasta`: Real-time growing list of representatives
- `unclustered_reads.fasta`: Reads that didn't cluster
- `clustering_summary.txt`: Comprehensive summary report

### Optional Output Files
- `cluster_consensus.fasta`: Consensus sequences (if `-p` used)
- `cluster_consensus_deduplicated.fasta`: Deduplicated consensus (if `--remove-redundant-consensus` used)
- `cluster_longest_reads.fasta`: Longest read from each cluster (if `--output-longest-reads` used)
- `cluster_*_consensus.fasta`: Individual consensus files for each cluster

## Configuration

### Parameter Tuning in `params.py`

The tool uses a centralized configuration file (`params.py`) for fine-tuning:

#### VSEARCH Parameters
```python
VSEARCH_PARAMS = {
    'clustering_algorithm': 'cluster_size',  # cluster_fast, cluster_size, cluster_unoise
    'identity_threshold': 0.9,               # 0.80-0.95
    'query_coverage': 0.9,                   # 0.70-0.90
    'target_coverage': 0.85,                 # 0.70-0.90
}
```

#### Redundant Removal Parameters
```python
REDUNDANT_REMOVAL_PARAMS = {
    'enabled': False,                        # Enable redundant removal
    'identity_threshold': 0.95,              # 95% identity threshold
    'query_coverage': 0.80,                  # 80% query coverage
    'target_coverage': 0.80,                 # 80% target coverage
}
```

## Notes on Algorithms

- Clustering is performed by VSEARCH (fixed).
- Consensus polishing, when enabled with `-p`, uses RACON + minimap2.

## Performance Tips

1. **For Large Datasets**: VSEARCH clustering is the default.
2. **For High Quality**: Enable RACON polishing with `-p`.
3. **For Memory Efficiency**: Set `--max-clusters` and `--max-read-length`
4. **For Debugging**: Use `--verbose --save-intermediate`

## Troubleshooting

### Common Issues

1. **VSEARCH not found**: Install VSEARCH and ensure it's in your PATH
2. **RACON/minimap2 errors**: Check that external tools are properly installed
3. **Memory issues**: Reduce `--max-clusters` or `--max-read-length`
4. **Slow performance**: Use VSEARCH clustering (`-f`) for large datasets

### Debug Mode

```bash
python cluster_polish.py input.fasta output_dir --debug --verbose --save-intermediate
```

## Citation

If you use this tool in your research, please cite:

```
Cluster Polish Tool - A modular read clustering pipeline with consensus polishing
Author: Assistant
Date: 2024
```

## License

This project is open source. Please check the license file for details.

## Contributing

Contributions are welcome! Please feel free to submit issues, feature requests, or pull requests.

## Support

For questions or issues, please:
1. Check the troubleshooting section
2. Review the example commands
3. Check that all dependencies are properly installed
4. Use debug mode for detailed error information
