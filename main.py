#!/usr/bin/env python3
"""Protein Sequence Homology Mapping Pipeline — main entry point.

Automates cross-referencing of a target species proteome against PDB
homologous structures within a parent taxonomic group.

Usage:
    python main.py --target 10253 --parent 10242 --keywords antibody
"""

from __future__ import annotations

import argparse
import sys
import time

import config
from taxonomy import verify_lineage
from uniprot import fetch_uniprot_proteins
from pdb_search import search_pdb_structures, download_all_pdb_fasta
from alignment import run_alignment
from clustering import run_mmseqs_clustering_report
from intermediate import (
    pdb_fasta_dir,
    save_pdb_chain_intermediates,
    save_pdb_search_intermediates,
    save_uniprot_intermediates,
)
from output import write_csv
from utils import log, probe_endpoints


def _write_empty_outputs_and_exit(
    target_taxid: int,
    parent_taxid: int,
    threshold: float,
    cluster_threshold: float,
    output_dir: str,
    target_name: str,
    target_proteins: list[dict[str, object]],
) -> None:
    """Write empty report files for a run that has no downstream candidates."""
    mapping_path, aln_path, log_path = write_csv([], [], target_taxid, parent_taxid, threshold, output_dir, target_name)
    cluster_path = run_mmseqs_clustering_report(
        [],
        target_proteins,
        target_taxid,
        parent_taxid,
        output_dir,
        cluster_threshold,
    )
    with open(log_path, "a", encoding="utf-8") as handle:
        handle.write(f"\nMMseqs2 Cluster Threshold: {cluster_threshold}\n")
        handle.write(f"Cluster Mapping CSV: {cluster_path}\n")
    log(f"  Empty mapping CSV: {mapping_path}")
    log(f"  Empty alignment detail CSV: {aln_path}")
    log(f"  Empty cluster mapping CSV: {cluster_path}")
    log(f"  Run log: {log_path}")
    sys.exit(0)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Protein sequence homology mapping pipeline: target proteome x PDB homologous structures",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--target", type=int, required=True, help="Target taxonomy ID")
    parser.add_argument("--parent", type=int, required=True, help="Parent taxonomy ID")
    parser.add_argument("--keywords", type=str, nargs="+", required=True, help="Space-separated PDB search keywords")
    parser.add_argument("--threshold", type=float, default=config.DEFAULT_THRESHOLD, help=f"Identity threshold %% (default: {config.DEFAULT_THRESHOLD})")
    parser.add_argument("--cluster-threshold", type=float, default=config.MMSEQS_CLUSTER_THRESHOLD, help=f"MMseqs2 clustering sequence identity threshold (default: {config.MMSEQS_CLUSTER_THRESHOLD})")
    parser.add_argument("--output", type=str, default=config.DEFAULT_OUTPUT_DIR, help=f"Output directory (default: {config.DEFAULT_OUTPUT_DIR})")
    args = parser.parse_args()

    target_taxid: int = args.target
    parent_taxid: int = args.parent
    keywords: list[str] = args.keywords
    threshold: float = args.threshold
    cluster_threshold: float = args.cluster_threshold
    output_dir: str = args.output

    start_time = time.time()
    log("=" * 60)
    log("Protein sequence homology mapping pipeline started")
    log(f"  Target Taxonomy ID: {target_taxid}")
    log(f"  Parent Taxonomy ID: {parent_taxid}")
    log(f"  Keywords: {keywords}")
    log(f"  Identity Threshold: {threshold}%")
    log(f"  MMseqs2 Cluster Threshold: {cluster_threshold}")
    log(f"  Output Directory: {output_dir}")
    log("=" * 60)

    # Probe all endpoint connectivity (cached for rest of pipeline)
    probe_endpoints()

    # 1. Verify taxonomy lineage
    log("\n[Step 1/6] Verify NCBI Taxonomy lineage")
    try:
        lineage_info = verify_lineage(target_taxid, parent_taxid)
        target_name = lineage_info["target_name"]
    except Exception as exc:
        log(f"  Lineage verification failed: {exc}", "ERROR")
        sys.exit(1)

    # 2. Download target proteome from UniProt
    log("\n[Step 2/6] Download target proteome from UniProt")
    try:
        target_proteins = fetch_uniprot_proteins(target_taxid)
        save_uniprot_intermediates(target_proteins, output_dir, target_taxid, parent_taxid)
    except Exception as exc:
        log(f"  UniProt download failed: {exc}", "ERROR")
        sys.exit(1)

    if not target_proteins:
        log("  UniProt returned no proteins; exiting", "WARNING")
        _write_empty_outputs_and_exit(
            target_taxid,
            parent_taxid,
            threshold,
            cluster_threshold,
            output_dir,
            target_name,
            [],
        )

    # 3. Search PDB structures
    log("\n[Step 3/6] Search RCSB PDB structures")
    try:
        pdb_ids = search_pdb_structures(parent_taxid, keywords)
        save_pdb_search_intermediates(pdb_ids, keywords, output_dir, target_taxid, parent_taxid)
    except Exception as exc:
        log(f"  PDB search failed: {exc}", "ERROR")
        sys.exit(1)

    if not pdb_ids:
        log("  No matching PDB structures found; exiting", "WARNING")
        _write_empty_outputs_and_exit(
            target_taxid,
            parent_taxid,
            threshold,
            cluster_threshold,
            output_dir,
            target_name,
            target_proteins,
        )

    # 4. Download PDB FASTA
    log("\n[Step 4/6] Download PDB FASTA sequences")
    try:
        raw_pdb_fasta_dir = pdb_fasta_dir(output_dir, target_taxid, parent_taxid)
        pdb_chains = download_all_pdb_fasta(pdb_ids, raw_pdb_fasta_dir)
        save_pdb_chain_intermediates(pdb_chains, output_dir, target_taxid, parent_taxid)
    except Exception as exc:
        log(f"  PDB FASTA download failed: {exc}", "ERROR")
        sys.exit(1)

    if not pdb_chains:
        log("  No PDB chain sequences were parsed; exiting", "WARNING")
        _write_empty_outputs_and_exit(
            target_taxid,
            parent_taxid,
            threshold,
            cluster_threshold,
            output_dir,
            target_name,
            target_proteins,
        )

    # 5. Sequence alignment
    log("\n[Step 5/6] Sequence alignment")
    aggregated, all_alignments = run_alignment(target_proteins, pdb_chains, threshold)

    # Write standard output files
    log("\n[Output] Write mapping and alignment CSV files")
    mapping_path, aln_path, log_path = write_csv(
        aggregated, all_alignments, target_taxid, parent_taxid, threshold, output_dir, target_name
    )

    # 6. Cluster mapped UniProt sequences
    log("\n[Step 6/6] Cluster mapped UniProt sequences with MMseqs2")
    try:
        cluster_path = run_mmseqs_clustering_report(
            aggregated,
            target_proteins,
            target_taxid,
            parent_taxid,
            output_dir,
            cluster_threshold,
        )
    except Exception as exc:
        log(f"  MMseqs2 clustering failed: {exc}", "ERROR")
        sys.exit(1)

    with open(log_path, "a", encoding="utf-8") as handle:
        handle.write(f"\nMMseqs2 Cluster Threshold: {cluster_threshold}\n")
        handle.write(f"Cluster Mapping CSV: {cluster_path}\n")

    elapsed = time.time() - start_time
    log("=" * 60)
    log(f"Pipeline completed in {elapsed:.1f}s")
    log(f"  Mapping CSV: {mapping_path}")
    log(f"  Alignment detail CSV: {aln_path}")
    log(f"  Cluster mapping CSV: {cluster_path}")
    log(f"  Run log: {log_path}")
    log(f"  Summary: {len(target_proteins)} target proteins, {len(pdb_ids)} PDB entries, {len(aggregated)} proteins with matches")
    log("=" * 60)


if __name__ == "__main__":
    main()
