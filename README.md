# Multiplex Primer Blast (MPB)

## Overview

### The Gap in Existing Tools

Tools like Primer3 and NCBI Primer-BLAST are excellent at designing and validating individual primer pairs, but they treat each pair in isolation. When you combine multiple pairs into a multiplex reaction, these tools have no concept of the pool as a whole — a forward primer from one pair is never checked against the reverse primers of all other pairs.

### Why This Matters at Scale

In a multiplex reaction, every primer can potentially interact with every other primer in the pool. With 10 primer pairs (20 primers), there are **190 possible primer combinations**. With 20 pairs, that number rises to **780**. No one reviews these manually, and no standard tool checks them automatically. Primers that each pass specificity checks individually can still produce unexpected bands, consume reagents, reduce sensitivity, or — in diagnostic settings — generate false positives. These problems only surface in the lab, after the primers are ordered.

### What Multiplex Primer Blast Does

**Multiplex Primer Blast (MPB)** fills this gap. Given a reference genome and a pool of primers, it aligns every primer against the reference simultaneously, then identifies every pair of binding sites that could produce a PCR product (a forward-strand hit and a reverse-strand hit within an amplicon-size window). Each candidate is filtered by an accurate thermodynamic Tm calculation, and every surviving amplicon is reported.

Because MPB makes no assumptions about which primers were designed as pairs, it doubles as a general-purpose specificity checker — you can use it for a single primer pair the same way you would use NCBI Primer-BLAST.

### Key Features
* **Pool-wide cross-amplicon detection** — every primer is checked against every other, not just within designed pairs.
* **Dual alignment backends** — NCBI BLAST+ (`blastn-short`) and MUMmer (`nucmer`). See [Choosing a Backend](#choosing-a-backend) for the tradeoff.
* **Thermodynamic filtering** — uses `primer3-py` to calculate accurate melting temperatures, filtering by Tm rather than by sequence identity alone.
* **Automatic BLAST DB management** — generates BLAST databases from raw reference FASTAs on first run.
* **Single-pair use** — works equally well for one primer pair as for a large multiplex pool.

---

## How It Works

1. **Align** every primer against the reference genome (BLAST or MUMmer).
2. **Classify** each hit as `forward` or `reverse` based on which strand it bound. *Note: this strand is derived from the alignment, not from the primer's intended role in your design — a primer you designed as a "reverse primer" will appear under `forward_primer` in any rows where it happened to bind the forward strand.*
3. **Calculate Tm** for each hit using `primer3-py` against the flanking reference sequence, and discard hits below `--tm-threshold`.
4. **Pair** every surviving forward hit with every reverse hit on the same chromosome where the distance falls within `[--min-amplicon, --max-amplicon]`.
5. **Report** every resulting amplicon to CSV.

### Reading the Output

MPB does not know which primers you designed as pairs, so it reports **every possible amplicon — including your intended on-target products**. This is deliberate, and it is what lets the same tool serve as both a multiplex cross-reactivity checker and a general-purpose primer specificity checker. In practice, your intended products are easy to recognise (expected size, expected chromosome, both primers from the same designed pair); everything else is a candidate off-target worth reviewing.

---

## Prerequisites

- **Python 3.8+**
- **primer3-py** — thermodynamic calculations
- **pandas** — parsing alignment results
- **BLAST+** *(recommended)* — `blastn` and `makeblastdb` on PATH
- **MUMmer4** *(alternative)* — `nucmer` and `show-coords` on PATH

### Installation via pip

`pip` only handles the Python dependencies. BLAST+ and MUMmer are external C binaries and must be installed separately so that `blastn`, `makeblastdb`, and `nucmer` are available on `PATH`.

**1. Install the alignment binaries.**

| Platform                 | BLAST+ and MUMmer install commands                                                       |
| ------------------------ | ---------------------------------------------------------------------------------------- |
| Ubuntu / Debian          | `sudo apt install ncbi-blast+ mummer`                                                    |
| Fedora / RHEL            | `sudo dnf install ncbi-blast+ mummer`                                                    |
| macOS (Homebrew)         | `brew install blast mummer`                                                              |
| Other / manual           | NCBI BLAST+: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ — MUMmer4: https://github.com/mummer4/mummer/releases |

You only need the backend you plan to use; install both if you want to switch with `--tool`.

**2. Install MPB and its Python dependencies.**

```bash
git clone https://github.com/uguraka/multiplex-primer-blast.git
cd multiplex-primer-blast

python -m venv mpb_env
source mpb_env/bin/activate   # on Windows: mpb_env\Scripts\activate

pip install -r requirements.txt
```

**3. Verify the binaries are on PATH:**

```bash
blastn -version
nucmer --version
```

### Installation via conda (recommended)

Installs MPB and all alignment binaries into a single environment.

```bash
git clone https://github.com/uguraka/multiplex-primer-blast.git
cd multiplex-primer-blast

conda create -n mpb_env -c bioconda python pandas primer3-py blast mummer
conda activate mpb_env
```

---

## Input Format

MPB takes two FASTA files:

- **`--ref`** — the reference genome (one or more sequences in a single FASTA).
- **`--primers`** — a **single multiFASTA containing every primer in the pool**, with one record per primer.

No naming convention is required, and there is no need to mark which primers are forward and which are reverse — MPB determines strand from the alignment. The names you choose appear verbatim in the output, so pick something meaningful.

Example primer FASTA:

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

---

## Usage

```bash
python MPB.py --ref <reference.fasta> --primers <primers.fasta> [options]
```

### Required Arguments
- `--ref` — path to the reference genome FASTA (e.g. `hg38_primary.fa`).
- `--primers` — path to the primer multiFASTA.

### Optional Arguments
- `--tool {blast,mummer}` — alignment backend. Default: `blast`.
- `--prefix` — prefix for output files. Default: `primer_alignment`.
- `--tm-threshold` — minimum melting temperature (°C) for a binding site to be considered valid. Default: `30.0`.
- `--min-amplicon` — minimum amplicon size in bp. Default: `50`.
- `--max-amplicon` — maximum amplicon size in bp. Default: `1000`.
- `--skip-alignment` — reuse an existing alignment file (useful when re-running with different Tm or amplicon-size thresholds without re-aligning).

### Reusing an Existing BLAST Database

The first time MPB runs with the BLAST backend, it builds a BLAST database from your reference FASTA and writes the `.nhr` / `.nin` / `.nsq` files alongside it (e.g. `hg38.fa.nhr`). On subsequent runs against the same FASTA, MPB detects these files and reuses the existing database instead of rebuilding it.

If you already have a pre-built BLAST database, place its `.nhr` / `.nin` / `.nsq` files next to the reference FASTA (or symlink them) using the FASTA filename as the prefix, and MPB will pick them up automatically. The FASTA itself is still required — MPB reads sequence from it during the Tm-calculation step.

### Choosing a Backend

|                              | BLAST (`blastn-short`)                                    | MUMmer (`nucmer`)                                  |
| ---------------------------- | --------------------------------------------------------- | -------------------------------------------------- |
| Recommended for              | most workflows                                            | very large references, when DB build is impractical |
| Sensitivity for short primers | tuned for short queries via `blastn-short`               | requires explicit short-match parameters           |
| Database setup               | first-run BLAST DB build, reused thereafter               | none — operates directly on the reference FASTA   |
| Output format                | tabular, easy to parse                                    | nucmer delta format                                |

BLAST is recommended for most users because `blastn-short` is purpose-built for short queries in the 18–25 bp range MPB targets, and it produces more sensitive hits for primer-length sequences. MUMmer is a reasonable alternative when the up-front cost of building a BLAST database is prohibitive.

---

## Output Format

MPB writes a CSV named `<prefix>_amplicon_results.csv` listing every candidate amplicon.

| column           | meaning                                                                 |
| ---------------- | ----------------------------------------------------------------------- |
| `chromosome`     | reference sequence on which the amplicon forms                          |
| `forward_primer` | name of the primer binding the forward strand (from `--primers`)        |
| `reverse_primer` | name of the primer binding the reverse strand (from `--primers`)        |
| `start`          | 1-based start position of the amplicon on the reference                 |
| `end`            | 1-based end position of the amplicon on the reference                   |
| `size`           | amplicon length in bp                                                   |
| `forward_tm`     | Tm (°C) of the forward primer at its binding site                       |
| `reverse_tm`     | Tm (°C) of the reverse primer at its binding site                       |
| `avg_tm`         | mean of `forward_tm` and `reverse_tm`                                   |

**Example output:**

| chromosome | forward_primer | reverse_primer | start  | end    | size | forward_tm | reverse_tm | avg_tm |
| ---------- | -------------- | -------------- | ------ | ------ | ---- | ---------- | ---------- | ------ |
| chr1       | Target1_F      | Target4_R      | 150230 | 150680 | 450  | 58.2       | 56.4       | 57.3   |
| chr5       | Target2_F      | Target2_R      | 884010 | 884110 | 100  | 62.1       | 61.8       | 61.95  |
| chrX       | Target9_F      | Target3_R      | 45200  | 45950  | 750  | 45.5       | 48.0       | 46.75  |

Row 1 is a cross-pair off-target. Row 2 is an intended on-target amplicon for the Target2 pair. Row 3 is another cross-pair off-target. MPB does not distinguish between these — interpretation is up to you.

---

## Limitations & Scope

- **Memory** — the full reference genome is loaded into a dictionary in memory. Multi-gigabyte references (e.g. human `hg38`) require sufficient RAM.
- **Plex size** — the number of primer-pair combinations grows quadratically with the number of primers (~N²/2). For pools above ~100 primers, expect longer runtimes during the thermodynamic-filtering step.
- **In-solution interactions are out of scope** — MPB only detects amplicons that form on the reference template. It does not check for primer dimers, hairpins, or other primer-primer interactions in solution. Use Primer3 or a dedicated dimer-checker for those.
- **Single reference per run** — MPB checks against one reference at a time. For panels intended to work across multiple genomes (e.g. host plus pathogen), run MPB separately against each reference.
