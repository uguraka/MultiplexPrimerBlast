# Multiplex Primer Blast (MPB)

## Overview

Designing primers for multiplex PCR is a complex challenge. While tools exist (like Primer3 or NCBI Primer-BLAST) to design highly specific single primer pairs, combining multiple pairs drastically increases the risk of **cross-reactivity** against the reference genome. 

When multiplexing, a forward primer from one pair might inadvertently interact with a reverse primer from an entirely different pair, forming **unwanted amplicons**. These unintended products consume PCR reagents, reduce amplification efficiency, and can produce spurious bands or false positives.

**Multiplex Primer Blast (MPB)** provides a command-line pipeline for analyzing the specificity of multiplex PCR primers against a full reference genome. It identifies potential unwanted amplicons by simulating the binding of an entire pool of primers and calculating the true thermodynamic melting temperature (Tm) of potential off-target binding sites. 

By exhaustively pairing valid forward and reverse binding sites within expected size constraints, MPB rapidly identifies unwanted amplicons that could ruin your multiplex reaction, allowing you to filter out problematic primers *before* ordering them.

### Key Features
* **Dual Alignment Backends:** Supports both NCBI BLAST+ (`blastn-short`) and MUMmer (`nucmer`) for rapid genomic alignment.
* **Automatic DB Management:** Automatically generates BLAST databases from raw reference FASTAs if they do not already exist.
* **Thermodynamic Filtering:** Uses `primer3-py` to calculate accurate binding melting temperatures, filtering out weak off-target alignments.

---

## Prerequisites

- **Python 3.8+**
- **primer3-py**: For thermodynamic calculations.
- **pandas**: For parsing and analyzing alignment results.
- **BLAST+**: (Default) Requires NCBI BLAST+ (`makeblastdb` and `blastn`) installed and in your system PATH.
- **MUMmer**: (Optional) Requires MUMmer4 (`nucmer`) installed and in your system PATH.

### Installation
```bash
# Clone the repository
git clone https://github.com/uguraka/MultiplexPrimerBlast.git
cd MultiplexPrimerBlast

# Create a virtual environment named 'mpb_env'
python -m venv mpb_env

# Activate the virtual environment
# On Linux/macOS:
source mpb_env/bin/activate
# On Windows:
# mpb_env\Scripts\activate

# Install the required Python packages
pip install -r requirements.txt

```

### Installing via conda

The most reliable way to install MPB and its dependencies (including the BLAST binaries) is via a Conda environment.

```bash
# Clone the repository
git clone https://github.com/uguraka/MultiplexPrimerBlast.git
cd MultiplexPrimerBlast

# Create and activate the conda environment
conda create -n mpb_env -c bioconda python pandas primer3-py blast
conda activate mpb_env
```

---



## Usage

### Command Line Interface

```bash
python MPB.py --ref <path_to_reference.fasta> --primers <path_to_primers.fasta> [options]
```

#### Required Arguments
- `--ref`: Path to the reference genome FASTA file (e.g., `hg38_primary.fa`).
- `--primers`: Path to the multiplex primer FASTA file.

#### Optional Arguments
- `--tool`: Alignment tool to use (`blast` or `mummer`). Default is `blast`.
- `--prefix`: Prefix for the output files. Default is `primer_alignment`.
- `--tm-threshold`: Minimum melting temperature (Tm) in °C for a binding site to be considered valid. Default is `30.0`.
- `--max-amplicon`: Maximum expected off-target amplicon size (in base pairs). Default is `1000`.
- `--min-amplicon`: Minimum expected off-target amplicon size (in base pairs). Default is `50`.
- `--skip-alignment`: Flag to skip the alignment step and use existing output files (useful for re-calculating Tm thresholds without re-running BLAST).

---

## Output Format

The script outputs a comprehensive CSV file named `<prefix>_amplicon_results.csv` detailing every potential unwanted amplicon discovered. 

**Example Output Data:**

| chromosome | forward_primer | reverse_primer | start | end | size | forward_tm | reverse_tm | avg_tm | priority |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| chr1 | Target1_F | Target4_R | 150230 | 150680 | 450 | 58.2 | 56.4 | 57.3 | |
| chr5 | Target2_F | Target2_R | 884010 | 884110 | 100 | 62.1 | 61.8 | 61.95 | |
| chrX | Target9_F | Target3_R | 45200 | 45950 | 750 | 45.5 | 48.0 | 46.75 | |


---

## Limitations & Performance

* **Memory Usage:** Analyzing large multi-gigabyte reference genomes (like human `hg38`) requires sufficient RAM to load sequence dictionaries and process large alignment tabular files. 
* **Plex-Size:** As the number of primers in the pool increases, the number of potential cross-reactive combinations increases exponentially. For highly complex pools (>100 plex), expect longer processing times during the thermodynamic calculation phase.

