"""Configuration constants for the protein homology mapping pipeline."""

from __future__ import annotations

# ── Output ─────────────────────────────────────────────────────────────────
DEFAULT_OUTPUT_DIR = "./results/"

# ── Network ────────────────────────────────────────────────────────────────
TIMEOUT = 60  # seconds for regular requests
FASTA_TIMEOUT = 30  # seconds for individual PDB FASTA downloads
PROBE_TIMEOUT = 15  # seconds for connectivity probe (per endpoint)
MAX_RETRIES = 3  # exponential backoff: 2s, 4s, 8s
RETRY_BACKOFF_BASE = 2.0  # seconds

# ── UniProt ────────────────────────────────────────────────────────────────
UNIPROT_PAGE_SIZE = 500

# ── RCSB PDB ───────────────────────────────────────────────────────────────
PDB_PAGE_SIZE = 1000
PDB_CONCURRENT_WORKERS = 5

# ── Alignment ──────────────────────────────────────────────────────────────
DEFAULT_THRESHOLD = 30.0  # percent identity
LENGTH_RATIO_MAX = 5.0  # skip if length ratio > this
KMER_SIZE = 3
KMER_JACCARD_THRESHOLD = 0.1
KMER_PRE_FILTER_TRIGGER = 10000  # total alignments threshold to enable k-mer pre-filter

# ── Intermediate outputs ───────────────────────────────────────────────────
SAVE_INTERMEDIATES = True
INTERMEDIATE_DIR_SUFFIX = "_intermediates"

# ── MMseqs2 clustering ─────────────────────────────────────────────────────
MMSEQS_BINARY = "mmseqs"
MMSEQS_CLUSTER_THRESHOLD = 0.9
MMSEQS_THREADS = 1
