# Installation

## Requirements

- Python 3.8+
- `primer3-py` and `pandas` (installed via pip)
- **BLAST+** (`blastn`, `makeblastdb`) or **MUMmer4** (`nucmer`, `show-coords`) — at least one required, depending on which backend you use

## Option 1 — conda (recommended)

Installs MPB and all alignment binaries into a single environment.

```bash
git clone https://github.com/uguraka/multiplex-primer-blast.git
cd multiplex-primer-blast

conda create -n mpb_env -c bioconda python pandas primer3-py blast mummer
conda activate mpb_env
```

## Option 2 — pip + manual binaries

**1. Install the alignment binaries.**

| Platform | BLAST+ | MUMmer4 |
|---|---|---|
| Ubuntu / Debian | `sudo apt install ncbi-blast+` | `sudo apt install mummer` |
| Fedora / RHEL | `sudo dnf install ncbi-blast+` | `sudo dnf install mummer` |
| macOS (Homebrew) | `brew install blast` | `brew install mummer` |
| Manual | [NCBI BLAST+ downloads](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) | [MUMmer4 releases](https://github.com/mummer4/mummer/releases) |

You only need the backend you plan to use.

**2. Clone and install Python dependencies.**

```bash
git clone https://github.com/uguraka/multiplex-primer-blast.git
cd multiplex-primer-blast

python -m venv mpb_env
source mpb_env/bin/activate   # Windows: mpb_env\Scripts\activate

pip install -r requirements.txt
```

**3. Verify the binaries are on PATH.**

```bash
blastn -version
nucmer --version
```

## Input format

MPB takes two FASTA files:

- **`--ref`** — the reference genome (one or more sequences in a single FASTA)
- **`--primers`** — a multiFASTA containing every primer in the pool, one record per primer

No naming convention is required. MPB determines strand from the alignment, not from primer names. Choose names that will be meaningful in the output CSV.

```
>Target1_F
ACGTACGTACGTACGTACGT
>Target1_R
TGCATGCATGCATGCATGCA
>Target2_F
GCATGCATGCATGCATGCAT
>Target2_R
CGTACGTACGTACGTACGTA
```
