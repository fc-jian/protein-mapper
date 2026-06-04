# Protein Sequence Homology Mapping Pipeline

[Chinese](README_CN.md)

Protein Mapper is a Python command-line pipeline for cross-species protein
sequence homology mapping. Given a target species, a parent taxonomic group, and
PDB search keywords, it fetches the target proteome from UniProt, searches RCSB
PDB structures under the parent lineage, downloads PDB FASTA chain sequences,
aligns sequences, and writes CSV mapping outputs.

## Workflow

```text
Input: target_taxid, parent_taxid, keywords[]
  |
  |-- [1] NCBI Datasets v2: verify parent lineage contains target
  |-- [2] UniProt REST API: fetch the target proteome
  |-- [3] RCSB PDB Search API v2: search parent-lineage structures by keywords
  |-- [4] RCSB PDB FASTA: download and parse PDB chain sequences
  |-- [5] PairwiseAligner + k-mer pre-screen: align sequences
  `-- [6] MMseqs2: cluster mapped UniProt sequences and write reports
```

## Requirements

- Python 3.10+
- `uv`
- MMseqs2 (`mmseqs` on `PATH`)
- Network access to NCBI, UniProt, and RCSB PDB public APIs

## Setup

This repository uses `uv` to manage the local environment in `./.venv`.

```bash
cd /path/to/protein-mapper
uv venv .venv --python 3.10
source .venv/bin/activate
uv pip install requests biopython pandas
```

Install MMseqs2 separately with your system or bioinformatics package manager,
then confirm that `mmseqs` is available on `PATH`.

## Run

```bash
source .venv/bin/activate
python main.py \
    --target 10253 \
    --parent 10242 \
    --keywords antibody Fab scFv nanobody VHH \
    --threshold 30.0 \
    --cluster-threshold 0.9 \
    --output ./results/
```

## CLI Arguments

| Argument | Required | Default | Description |
| --- | :---: | --- | --- |
| `--target` | yes | - | Target NCBI taxonomy ID |
| `--parent` | yes | - | Parent NCBI taxonomy ID used for lineage filtering |
| `--keywords` | yes | - | Space-separated PDB search keywords |
| `--threshold` | no | `30.0` | Minimum identity percentage |
| `--cluster-threshold` | no | `0.9` | MMseqs2 sequence identity threshold for clustering mapped UniProt proteins |
| `--output` | no | `./results/` | Output directory |

## Outputs

Each run writes four report files:

| File | Description |
| --- | --- |
| `{target}_{parent}_mapping.csv` | Main mapping table, one row per matched target protein |
| `{target}_{parent}_all_alignments.csv` | All alignments passing the identity threshold |
| `{target}_{parent}_cluster_mapping.csv` | Cluster-level mapping table after MMseqs2 clustering of mapped UniProt sequences |
| `{target}_{parent}_log.txt` | Summary log for the run |

`results/` is treated as generated output and is not tracked by Git.

Intermediate sequence and metadata files are saved by default:

```text
{output}/{target}_{parent}_intermediates/
├── uniprot/
│   ├── {target}_uniprot.fasta
│   └── {target}_uniprot_proteins.csv
├── pdb/
│   ├── {target}_{parent}_pdb_ids.txt
│   ├── {target}_{parent}_pdb_query.json
│   ├── {target}_{parent}_pdb_chains.csv
│   └── fasta/{pdb_id}.fasta
└── mmseqs/
    └── {target}_{parent}_mapped_uniprot.fasta
```

The cluster report groups final mapped UniProt sequences at the configured
MMseqs2 identity threshold. For each cluster, all mapping columns are
deduplicated and joined with semicolons.

### Mapping Columns

| Column | Example |
| --- | --- |
| `target_taxonomy_id` | `10253` |
| `parent_taxonomy_id` | `10242` |
| `target_protein_uniprot_id` | `Q9JFA1` |
| `target_protein_name` | `Cell surface-binding protein OPG105` |
| `target_organism` | `Vaccinia virus (strain Tian Tan)` |
| `support_pdb_id` | `4E9O;4ETQ;5USH;5USL;6B9J` |
| `support_pdb_chain_description` | PDB FASTA chain descriptions |
| `identity_pct` | `83.28` |
| `alignment_score` | `259.0` |

When multiple PDB entries match the same target protein, supporting PDB IDs and
chain descriptions are joined with semicolons. For each PDB ID, the best
identity chain is retained in the aggregated mapping.

## Example

```bash
python main.py --target 10253 --parent 10242 --keywords antibody --output ./results/
```

Example summary from a successful run:

```text
[Step 1/6] Verify NCBI Taxonomy lineage
[Step 2/6] Download target proteome from UniProt
[Step 3/6] Search RCSB PDB structures
[Step 4/6] Download PDB FASTA sequences
[Step 5/6] Sequence alignment
[Step 6/6] Cluster mapped UniProt sequences with MMseqs2

Pipeline completed.
```

## API Endpoints

| Step | API | Method | Endpoint |
| --- | --- | --- | --- |
| 1 | NCBI Datasets v2 | GET | `api.ncbi.nlm.nih.gov/datasets/v2/taxonomy/taxon/{taxid}` |
| 2 | UniProt REST | GET | `rest.uniprot.org/uniprotkb/search?query=taxonomy_id:{taxid}` |
| 3 | RCSB PDB Search v2 | POST | `search.rcsb.org/rcsbsearch/v2/query` |
| 4 | RCSB PDB FASTA | GET | `www.rcsb.org/fasta/entry/{pdb_id}/download` |

## Repository Layout

```text
protein-mapper/
├── main.py        # CLI entry point and pipeline orchestration
├── config.py      # Timeout, pagination, concurrency, and threshold constants
├── utils.py       # HTTP endpoint probing, retry handling, and logging
├── taxonomy.py    # NCBI lineage validation
├── uniprot.py     # UniProt proteome download and JSON parsing
├── pdb_search.py  # RCSB PDB search and FASTA download/parsing
├── alignment.py   # PairwiseAligner plus k-mer pre-screening
├── intermediate.py # Intermediate FASTA and metadata writers
├── clustering.py  # MMseqs2 clustering and cluster-level report generation
├── output.py      # CSV and summary-log writers
├── AGENTS.md      # Contributor and agent guidance
├── README.md      # English documentation
└── README_CN.md   # Chinese documentation
```

## Configuration

Common tuning constants live in `config.py`:

```python
TIMEOUT = 60
FASTA_TIMEOUT = 30
MAX_RETRIES = 3
SAVE_INTERMEDIATES = True
INTERMEDIATE_DIR_SUFFIX = "_intermediates"
UNIPROT_PAGE_SIZE = 500
PDB_PAGE_SIZE = 1000
PDB_CONCURRENT_WORKERS = 5
DEFAULT_THRESHOLD = 30.0
LENGTH_RATIO_MAX = 5.0
KMER_SIZE = 3
KMER_JACCARD_THRESHOLD = 0.1
KMER_PRE_FILTER_TRIGGER = 10000
MMSEQS_BINARY = "mmseqs"
MMSEQS_CLUSTER_THRESHOLD = 0.9
MMSEQS_THREADS = 1
```

## Dependencies

| Package | Purpose |
| --- | --- |
| `requests` | HTTP requests |
| `biopython` | Pairwise alignment and FASTA parsing |
| `pandas` | CSV output |
| MMseqs2 | Clustering mapped UniProt sequences |

## License

Internal research tool.
