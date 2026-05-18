# Usage

## Basic command

```bash
python MPB.py --ref <reference.fasta> --primers <primers.fasta> [options]
```

## Arguments

### Required

| Argument | Description |
|---|---|
| `--ref` | Path to the reference genome FASTA |
| `--primers` | Path to the primer multiFASTA |

### Optional

| Argument | Default | Description |
|---|---|---|
| `--tool {blast,mummer}` | `blast` | Alignment backend |
| `--prefix` | `primer_alignment` | Prefix for output files |
| `--tm-threshold` | `30.0` | Minimum Tm (°C) for a binding site to be reported |
| `--min-amplicon` | `50` | Minimum amplicon size in bp |
| `--max-amplicon` | `1000` | Maximum amplicon size in bp |
| `--skip-alignment` | — | Reuse an existing alignment file |

## Examples

Run with BLAST (default):
```bash
python MPB.py --ref hg38.fa --primers my_panel.fa
```

Run with MUMmer:
```bash
python MPB.py --ref hg38.fa --primers my_panel.fa --tool mummer
```

Tighten the Tm threshold and amplicon window:
```bash
python MPB.py --ref hg38.fa --primers my_panel.fa --tm-threshold 45 --min-amplicon 80 --max-amplicon 400
```

Re-run with different thresholds without re-aligning:
```bash
python MPB.py --ref hg38.fa --primers my_panel.fa --skip-alignment --tm-threshold 50
```

## Choosing a backend

| | BLAST (`blastn-short`) | MUMmer (`nucmer`) |
|---|---|---|
| Recommended for | most workflows | very large references where DB build is impractical |
| Short-primer sensitivity | tuned via `blastn-short` | requires explicit short-match parameters |
| Database setup | first-run build, reused thereafter | none — operates on the FASTA directly |

BLAST is recommended for most users. `blastn-short` is purpose-built for 18–25 bp queries.

## Output

MPB writes `<prefix>_amplicon_results.csv`:

| Column | Description |
|---|---|
| `chromosome` | Reference sequence on which the amplicon forms |
| `forward_primer` | Name of the primer binding the forward strand |
| `reverse_primer` | Name of the primer binding the reverse strand |
| `start` | 1-based start position of the amplicon |
| `end` | 1-based end position of the amplicon |
| `size` | Amplicon length in bp |
| `forward_tm` | Tm (°C) of the forward primer at its binding site |
| `reverse_tm` | Tm (°C) of the reverse primer at its binding site |
| `avg_tm` | Mean of `forward_tm` and `reverse_tm` |

MPB reports **every possible amplicon — including your intended on-target products**. It makes no assumptions about which primers were designed as pairs. Your intended products are recognisable by expected size, chromosome, and primer names; everything else is a candidate off-target.

## BLAST database reuse

On first run, MPB builds a BLAST database from your reference FASTA (`.nhr`, `.nin`, `.nsq` files alongside the reference). On subsequent runs against the same FASTA, the existing database is reused automatically. If you have a pre-built database, place its files next to the reference FASTA using the FASTA filename as the prefix.

## Limitations

- **Memory** — the full reference is loaded into memory. Multi-gigabyte references (e.g. hg38) require sufficient RAM.
- **Plex size** — runtime during the Tm-filtering step grows with the number of primer combinations (~N²/2). For pools above ~100 primers, expect longer runtimes.
- **In-solution interactions** — MPB only detects amplicons that form on the reference template. Primer dimers, hairpins, and other solution-phase interactions are out of scope; use Primer3 for those.
- **Single reference per run** — for panels intended to work across multiple genomes, run MPB separately against each reference.
