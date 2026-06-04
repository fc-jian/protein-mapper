---
name: mmseqs-cluster-report
description: Use when implementing, reviewing, or changing MMseqs2 clustering of mapped UniProt sequences, mapped FASTA generation, cluster TSV parsing, or cluster-level CSV reporting in the protein-mapper project.
---

# MMseqs2 Cluster Report

## Purpose

Cluster the final mapped UniProt sequences after sequence alignment and produce
a cluster-level CSV report. This step runs after the standard mapping and
all-alignments CSV files are written.

## Inputs

- `aggregated`: final mapping rows from `alignment.run_alignment`.
- `target_proteins`: normalized UniProt protein dictionaries.
- Default sequence identity threshold: `config.MMSEQS_CLUSTER_THRESHOLD` (`0.9`).

Only target proteins that passed final mapping are written to the MMseqs2 input
FASTA.

## MMseqs2 Command

Use `mmseqs easy-cluster`:

```text
mmseqs easy-cluster mapped_uniprot.fasta output_prefix tmp --min-seq-id 0.9 --threads N
```

Read the cluster membership from `{output_prefix}_cluster.tsv`, where each line
contains representative and member accessions.

## Output Rules

Write `{target}_{parent}_cluster_mapping.csv`.

Cluster report columns start with:

- `cluster_id`
- `cluster_representative`
- `cluster_members`
- `cluster_size`

Then include every standard mapping column. For every mapping column, merge all
cluster-member values by splitting existing semicolon-delimited values,
deduplicating while preserving first-seen order, and joining with semicolons.

## Failure Handling

If `mmseqs` is missing or exits nonzero, fail loudly with the MMseqs2 error.
Do not silently replace MMseqs2 clustering with a different algorithm.
