# How It Works

## Pipeline

MPB runs a five-step pipeline:

```
primers.fasta + reference.fasta
        │
        ▼
1. ALIGN — run BLAST or MUMmer to find all primer binding sites
        │
        ▼
2. CLASSIFY — label each hit as forward (+) or reverse (−) strand
        │
        ▼
3. FILTER — calculate Tm with primer3-py; drop hits below --tm-threshold
        │
        ▼
4. PAIR — find all (forward, reverse) pairs on the same chromosome
          within [--min-amplicon, --max-amplicon]
        │
        ▼
5. REPORT — write amplicons to CSV
```

## Alignment backends

### BLAST (`blastn-short`)

MPB calls `blastn` with the `blastn-short` task, which is tuned for short queries in the 18–25 bp range. On first run it builds a BLAST database from the reference FASTA using `makeblastdb`; the database files are reused on subsequent runs.

### MUMmer (`nucmer`)

MPB calls `nucmer` followed by `show-coords` to produce a tabular coordinates file. No database build is required, making MUMmer a practical choice for very large references.

Both backends produce the same internal `PrimerHit` records, so everything downstream is identical regardless of which backend was used.

## Thermodynamic filtering

For each alignment hit, MPB extracts the flanking reference sequence around the binding site and calls `primer3.bindings.calc_end_stability()` to compute a Tm. Hits below `--tm-threshold` (default 30 °C) are discarded. This step filters out low-affinity alignments that would never produce a product in practice — something sequence-identity cutoffs alone cannot do.

## Amplicon pairing algorithm

The pairing step needs to find every (forward hit, reverse hit) pair on the same chromosome where the distance falls within the amplicon size window. A naïve O(F × R) search would be too slow for large primer pools against whole genomes.

MPB uses an O(F + R) sliding-window approach:

1. Sort forward hits by start position.
2. Sort reverse hits by end position.
3. Advance a left pointer to drop forward hits whose start position is too far left to pair with the current reverse hit.
4. Emit all remaining forward hits that fall within the window.

This means the algorithm scales linearly with the number of hits rather than quadratically.

## Strand classification

Each alignment hit is classified based on which strand it bound — not based on the primer's name or intended role. A primer you designed as a "reverse primer" will appear under `forward_primer` in any rows where it happened to bind the forward strand of the reference. This is intentional: MPB checks all primer-primer combinations regardless of design intent.
