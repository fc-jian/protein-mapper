# AGENTS.md

## Scope

These instructions apply to the whole repository.

## Project Overview

This is a Python command-line pipeline for protein sequence homology mapping.
Given a target taxonomy ID, a parent taxonomy ID, and PDB search keywords, the
pipeline:

1. Verifies with NCBI Taxonomy that the parent taxon is in the target lineage.
2. Fetches the target proteome from UniProt as paginated JSON.
3. Searches RCSB PDB for structures under the parent taxonomy lineage with the
   requested keyword text query.
4. Downloads each matching PDB entry's FASTA records and parses chains.
5. Saves intermediate FASTA/metadata files by default.
6. Aligns target proteins against PDB chains and writes CSV outputs.
7. Clusters mapped UniProt sequences with MMseqs2 and writes a cluster report.

Primary tracked documentation lives in `README.md` and `README_CN.md`. The
`docs/` directory is local-only and ignored by Git.

## Local Skills

Reusable project knowledge is persisted under `.agents/skills/`:

- `.agents/skills/uniprot-proteome-query/SKILL.md`: UniProt taxonomy queries,
  pagination, and normalized target protein records.
- `.agents/skills/rcsb-pdb-query/SKILL.md`: RCSB Search API queries, taxonomy
  lineage filters, FASTA downloads, and chain parsing.
- `.agents/skills/protein-sequence-alignment/SKILL.md`: pairwise alignment,
  identity calculation, k-mer pre-screening, filtering, and aggregation.
- `.agents/skills/mmseqs-cluster-report/SKILL.md`: MMseqs2 clustering of mapped
  UniProt sequences and cluster-level CSV reporting.

When a task touches one of these areas, read the matching skill before editing.

## Runtime

- Python 3.10+.
- Runtime dependencies: `requests`, `biopython`, and `pandas`.
- External runtime dependency: MMseqs2 (`mmseqs` on `PATH`) for clustering the
  final mapped UniProt sequences.
- Use `uv` for the local environment in `./.venv`.
- There is no dependency manifest yet. If dependencies change, add or update a
  manifest and keep the README files in sync.

Setup:

```bash
uv venv .venv --python 3.10
source .venv/bin/activate
uv pip install requests biopython pandas
```

Run the documented smoke case:

```bash
python main.py --target 10253 --parent 10242 --keywords antibody --threshold 30.0 --output ./results/
```

Cheap local syntax check:

```bash
python -m py_compile main.py config.py utils.py taxonomy.py uniprot.py pdb_search.py alignment.py output.py
```

## Network And API Constraints

- The pipeline needs network access to NCBI, UniProt, and RCSB PDB public APIs.
- Use the shared `utils.http_get` and `utils.http_post` helpers for outbound
  calls so timeout, retry, and endpoint-probe behavior stays consistent.
- Keep NCBI lineage validation before the heavier UniProt and RCSB downloads.
- RCSB PDB FASTA is downloaded one PDB entry at a time and parsed with
  `Bio.SeqIO`; do not assume a batch FASTA endpoint exists.
- Keep timeout, retry, page-size, threshold, and concurrency knobs in
  `config.py`. Be conservative when changing `PDB_CONCURRENT_WORKERS`, because
  remote API rate behavior and local network throughput can dominate runtime.

Main endpoints used by the pipeline:

- NCBI Datasets v2 taxonomy: `GET /datasets/v2/taxonomy/taxon/{taxid}`
- UniProt REST search: `GET /uniprotkb/search?query=taxonomy_id:{taxid}`
- RCSB Search API v2: `POST /rcsbsearch/v2/query`
- RCSB entry FASTA: `GET /fasta/entry/{pdb_id}/download`

## Code Guidelines

- Keep the current small-module layout:
  - `main.py`: CLI parsing and orchestration.
  - `config.py`: constants and tuning knobs.
  - `utils.py`: HTTP endpoint probes, retry, and logging helpers.
  - `taxonomy.py`: NCBI lineage verification.
  - `uniprot.py`: UniProt proteome download and JSON parsing.
  - `pdb_search.py`: RCSB search plus concurrent FASTA download/parsing.
  - `alignment.py`: sequence alignment and aggregation.
  - `intermediate.py`: intermediate FASTA and metadata writers.
  - `clustering.py`: MMseqs2 clustering and cluster-level report generation.
  - `output.py`: CSV and run-log writing.
- Follow the existing style: type hints, `from __future__ import annotations`,
  module-level docstrings, and simple dictionaries for API records.
- Prefer structured JSON fields from UniProt and `Bio.SeqIO` for FASTA parsing.
  Avoid fragile string parsing unless the source format is documented and local.
- The current alignment implementation uses `Bio.Align.PairwiseAligner` in local
  mode, a length-ratio filter, and optional k-mer pre-screening for large
  all-vs-all jobs. Keep docs and comments aligned if changing this behavior.
- Preserve the CLI contract and output column names unless the user explicitly
  asks for a breaking change.

## Outputs And Data

The pipeline writes three files per run:

- `{target}_{parent}_mapping.csv`
- `{target}_{parent}_all_alignments.csv`
- `{target}_{parent}_cluster_mapping.csv`
- `{target}_{parent}_log.txt`

Intermediate FASTA and metadata files are saved by default under
`{output}/{target}_{parent}_intermediates/`.

`results/` and `*_intermediates/` are generated output and ignored by Git. Do
not force-add generated datasets unless the task specifically requires it.
`output.py` currently writes CSV files with `utf-8-sig`; preserve that unless
there is a concrete compatibility reason to change it.

## Validation

There is no formal test suite in this repository right now.

For small edits, run the local syntax check above. For changes to API clients,
FASTA parsing, alignment, or output schema, add focused fixture tests where
practical and run the documented smoke case when network access is available. If
a live API smoke run is skipped, state that explicitly in the handoff.

## Documentation

When changing CLI arguments, dependencies, endpoint behavior, output columns, or
alignment semantics, update both `README.md` and `README_CN.md` in the same
change.
