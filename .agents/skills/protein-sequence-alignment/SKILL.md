---
name: protein-sequence-alignment
description: Use when implementing, reviewing, or changing protein sequence alignment, identity calculation, k-mer pre-screening, length filtering, or alignment result aggregation in the protein-mapper project.
---

# Protein Sequence Alignment

## Purpose

Run all-vs-all local alignments between normalized target proteins and parsed
PDB chains, then keep alignments whose identity percentage meets the configured
threshold.

## Current Inputs

Target protein dictionaries are produced by the UniProt layer and must include:

- `accession`
- `entry_name`
- `protein_name`
- `sequence`
- `length`
- `organism`
- `taxon_id`

PDB chain dictionaries are produced by the RCSB layer and must include:

- `pdb_id`
- `chain_desc`
- `sequence`
- `length`

## Pre-Filters

Apply cheap filters before pairwise alignment:

- Length ratio: skip pairs where `max(lengths) / min(lengths)` exceeds
  `config.LENGTH_RATIO_MAX`.
- k-mer Jaccard: when total pair count exceeds
  `config.KMER_PRE_FILTER_TRIGGER`, compute `KMER_SIZE`-mers and skip pairs
  below `config.KMER_JACCARD_THRESHOLD`.

Keep these thresholds in `config.py`.

## Alignment

The current implementation uses `Bio.Align.PairwiseAligner` in local mode.
Identity is calculated from the best alignment as:

```text
matching aligned positions / aligned length * 100
```

Return `(identity_percent, alignment_score)`. If no alignment exists, return
`(0.0, 0.0)`.

When changing scoring, state the scoring scheme explicitly in code comments and
documentation. If a substitution matrix such as BLOSUM62 is intended, set it on
the aligner directly; do not rely on implicit defaults.

## Result Records

Each passing alignment should include:

- target accession, entry name, protein name, sequence, and length
- PDB ID, chain description, sequence, and length
- rounded `identity_pct`
- rounded `alignment_score`

Aggregation is one row per target protein:

- group by target accession
- keep the highest-identity chain per PDB ID
- join supporting PDB IDs and chain descriptions with semicolons
- sort final rows by descending identity

## Validation

Prefer deterministic unit tests for:

- `kmer_jaccard`
- length-ratio filtering behavior
- identity threshold inclusion/exclusion
- aggregation choosing the best chain per PDB

Avoid relying on live API calls for alignment tests; use small synthetic
sequence fixtures.

## Downstream Clustering

After final mapping aggregation, mapped UniProt sequences are clustered with
MMseqs2 at `config.MMSEQS_CLUSTER_THRESHOLD` (default `0.9`). Do not change the
mapping row schema without updating cluster report generation as well.
