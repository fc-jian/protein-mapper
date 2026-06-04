---
name: uniprot-proteome-query
description: Use when implementing, reviewing, or changing UniProt REST proteome queries by taxonomy ID, paginated result fetching, UniProt JSON parsing, or target protein record schemas in the protein-mapper project.
---

# UniProt Proteome Query

## Purpose

Use UniProtKB REST search to fetch all protein entries for one NCBI taxonomy ID
as structured JSON records. Keep this path JSON-first so downstream CSV output
has stable accessions, names, organisms, lengths, and sequences.

## Query Pattern

Use `GET https://rest.uniprot.org/uniprotkb/search` with:

- `query=taxonomy_id:{taxid}`
- `format=json`
- `size=500`

In this repository, call the shared HTTP helper instead of raw `requests`:

```python
resp = http_get(url)
data = resp.json()
```

## Pagination

UniProt returns the next page in the `Link` response header with
`rel="next"`. Continue until no next link remains.

Preferred loop shape:

```python
url = base_url
all_results = []
while url:
    resp = http_get(url)
    data = resp.json()
    all_results.extend(data.get("results", []))
    url = None
    for part in resp.headers.get("Link", "").split(","):
        if 'rel="next"' in part:
            url = part.split(";")[0].strip(" <>")
            break
```

## Target Protein Record

Normalize UniProt entries to dictionaries with these keys:

- `accession`: `primaryAccession`
- `entry_name`: `uniProtkbId`
- `protein_name`: recommended full name, falling back to first submitted name
- `sequence`: `sequence.value`
- `length`: `sequence.length`
- `organism`: `organism.scientificName`
- `taxon_id`: `organism.taxonId`

Avoid parsing FASTA headers for these fields when JSON is available.

## Edge Cases

- Missing recommended names are normal; fall back to submitted names, then `""`.
- Missing or empty `results` should return an empty list, not crash.
- Do not assume a fixed protein count for a taxonomy ID in tests; use small
  fixtures for parser tests and treat live API counts as smoke-test observations.
- Keep the output schema synchronized with `alignment.py`, `output.py`, and the
  README whenever normalized keys change.

## Validation

For parser changes, add fixture-level tests around `_parse_uniprot_results`.
For client changes, validate pagination with a mocked `Link` header before using
live API smoke runs.

## Intermediate Output

The pipeline saves normalized UniProt records by default under
`{output}/{target}_{parent}_intermediates/uniprot/`:

- `{target}_uniprot.fasta`
- `{target}_uniprot_proteins.csv`

Keep these files synchronized with the normalized record schema.
