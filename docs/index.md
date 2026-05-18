# Multiplex Primer Blast

**Validate PCR primer pools by detecting unintended cross-reactivity amplicons.**

Tools like Primer3 and NCBI Primer-BLAST are excellent at designing and validating individual primer pairs, but they treat each pair in isolation. When you combine multiple pairs into a multiplex reaction, a forward primer from one pair is never checked against the reverse primers of all other pairs.

**Multiplex Primer Blast (MPB)** fills this gap. Given a reference genome and a pool of primers, it aligns every primer against the reference simultaneously, then identifies every pair of binding sites that could produce a PCR product. Each candidate is filtered by an accurate thermodynamic Tm calculation, and every surviving amplicon is reported.

## Why it matters

In a multiplex reaction, every primer can potentially interact with every other primer in the pool. With 10 primer pairs (20 primers), there are **190 possible primer combinations**. With 20 pairs, that number rises to **780**. No one reviews these manually, and no standard tool checks them automatically.

Primers that each pass specificity checks individually can still produce unexpected bands, consume reagents, reduce sensitivity, or — in diagnostic settings — generate false positives. These problems only surface in the lab, after the primers are ordered.

## Key features

- **Pool-wide cross-amplicon detection** — every primer is checked against every other, not just within designed pairs
- **Dual alignment backends** — NCBI BLAST+ (`blastn-short`) and MUMmer (`nucmer`)
- **Thermodynamic filtering** — uses `primer3-py` to calculate accurate melting temperatures
- **Automatic BLAST DB management** — generates and reuses BLAST databases automatically
- **Single-pair use** — works equally well for one primer pair as for a large multiplex pool

## Quick start

```bash
git clone https://github.com/uguraka/multiplex-primer-blast.git
cd multiplex-primer-blast
pip install -r requirements.txt

python MPB.py --ref reference.fasta --primers primers.fasta
```

See [Installation](installation.md) for full setup including BLAST+ and MUMmer binaries.
