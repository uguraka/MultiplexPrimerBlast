# Multiplex Primer Blast (MPB)

## Overview

### The Gap in Existing Tools

Tools like Primer3 and NCBI Primer-BLAST are excellent at designing specific primer pairs — but they are built around a single pair at a time. When you combine multiple pairs into a multiplex reaction, these tools have no concept of the pool as a whole. A forward primer from one pair is never checked against the reverse primers of all other pairs. This is not a flaw — it is simply outside their scope.

### Why This Matters at Scale

In a multiplex reaction, every primer can potentially interact with every other primer in the pool. With 10 primer pairs (20 primers), there are **190 possible cross-pair combinations**. With 20 pairs, that number rises to **780**. No one reviews these manually, and no standard tool checks them automatically. The result: primers that each pass specificity checks individually can still produce unexpected bands, consume reagents, reduce sensitivity, or — in diagnostic settings — generate false positives. These problems only surface in the lab, after the primers are ordered and the reaction is running.

### What MPB Does

**Multiplex Primer Blast (MPB)** fills this gap. Given a reference genome and a pool of primers, it exhaustively checks every cross-pair combination by aligning all primers simultaneously and simulating their thermodynamic binding. It filters out weak off-target interactions using accurate Tm calculations, then reports every potential unwanted amplicon — before you order a single primer.

### Key Features
* **Dual Alignment Backends:** Supports both NCBI BLAST+ (`blastn-short`) and MUMmer (`nucmer`) for rapid genomic alignment.
* **Automatic DB Management:** Automatically generates BLAST databases from raw reference FASTAs if they do not already exist.
* **Thermodynamic Filtering:** Uses `primer3-py` to calculate accurate binding melting temperatures, filtering out weak off-target alignments.

---

## Prerequisites

- **Python 3.8+**
- **primer3-py**: For thermodynamic calculations.
- **pandas**: For parsing and analyzing alignment results.
- **BLAST+**: (Recommended) Requires NCBI BLAST+ (`makeblastdb` and `blastn`) installed and in your system PATH.
- **MUMmer**: (Alternative) Requires MUMmer4 (`nucmer`) installed and in your system PATH.

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
conda create -n mpb_env -c bioconda python pandas primer3-py blast mummer
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
- `--tool`: Alignment tool to use (`blast` or `mummer`). Default is `blast` (recommended).
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

