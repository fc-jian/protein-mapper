---
name: rcsb-pdb-query
description: Use when implementing, reviewing, or changing RCSB PDB Search API v2 queries, taxonomy lineage filters, keyword search behavior, PDB FASTA downloading, or PDB chain parsing in the protein-mapper project.
---

# RCSB PDB Query

## Purpose

Find PDB entries associated with a parent taxonomy lineage and parse each
entry's FASTA records into chain-level sequence dictionaries for alignment.

## Search API

Use `POST https://search.rcsb.org/rcsbsearch/v2/query` with
`return_type: "entry"` and paginated `request_options`.

Core taxonomy filter:

```json
{
  "type": "terminal",
  "service": "text",
  "parameters": {
    "attribute": "rcsb_entity_source_organism.taxonomy_lineage.id",
    "operator": "exact_match",
    "value": "10242"
  }
}
```

Keyword search should stay explicit about semantics. The current project uses
text search over `struct.title` and `rcsb_polymer_entity.pdbx_description`,
joined with an OR group. If strict all-keyword behavior is required, validate
RCSB operator semantics and prefer one condition per keyword.

Use `config.PDB_PAGE_SIZE` and continue paging until `start >= total_count` or
the API returns an empty batch.

## FASTA Download

Download one FASTA file per PDB ID:

```text
GET https://www.rcsb.org/fasta/entry/{PDB_ID}/download
```

There is no assumed batch FASTA endpoint in this project. Keep concurrency
bounded by `config.PDB_CONCURRENT_WORKERS`.

## Chain Record

Parse FASTA with `Bio.SeqIO` from an in-memory `StringIO`. Normalize non-empty
records to dictionaries with:

- `pdb_id`: entry ID from the search result
- `chain_desc`: full FASTA description line
- `sequence`: amino acid sequence
- `length`: sequence length

Do not split the chain description unless a caller needs a specific subfield.
The complete description is preserved in CSV output for traceability.

## Failure Handling

- Individual FASTA download failures should log a warning and skip that PDB.
- FASTA parser failures for one PDB should not stop the entire run.
- Empty search results and empty parsed chains should propagate as empty lists so
  `main.py` can write empty output files with headers.

## Validation

Use small FASTA fixtures to test parsing. For search changes, inspect the JSON
query body and use a live smoke run only when network access is available.

## Intermediate Output

The pipeline saves PDB search and FASTA intermediates by default under
`{output}/{target}_{parent}_intermediates/pdb/`:

- `{target}_{parent}_pdb_ids.txt`
- `{target}_{parent}_pdb_query.json`
- `{target}_{parent}_pdb_chains.csv`
- `fasta/{pdb_id}.fasta`

Raw FASTA files should be written after successful downloads and before parsing.
