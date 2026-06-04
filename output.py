"""CSV output writer for the protein homology mapping pipeline."""

from __future__ import annotations

import csv
import os
from typing import Any

import pandas as pd

from utils import log


def write_csv(
    aggregated: list[dict[str, Any]],
    all_alignments: list[dict[str, Any]],
    target_taxid: int,
    parent_taxid: int,
    thresholds: float,
    output_dir: str = ".",
    target_name: str = "",
) -> tuple[str, str, str]:
    """Write pipeline results to CSV files.

    Produces three files:
        1. {target}_{parent}_mapping.csv — aggregated summary (one row per target protein)
        2. {target}_{parent}_all_alignments.csv — all passing alignments
        3. {target}_{parent}_log.txt — brief summary log

    Args:
        aggregated: Aggregated results (one per target protein).
        all_alignments: All individual alignment records.
        target_taxid: Target taxonomy ID.
        parent_taxid: Parent taxonomy ID.
        thresholds: Identity threshold used.
        output_dir: Output directory path.
        target_name: Human-readable target organism name.

    Returns:
        Tuple of (mapping_path, alignments_path, log_path).
    """
    os.makedirs(output_dir, exist_ok=True)
    base = f"{target_taxid}_{parent_taxid}"
    mapping_path = os.path.join(output_dir, f"{base}_mapping.csv")
    alignments_path = os.path.join(output_dir, f"{base}_all_alignments.csv")
    log_path = os.path.join(output_dir, f"{base}_log.txt")

    # Fill parent_taxonomy_id
    for row in aggregated:
        row["parent_taxonomy_id"] = parent_taxid

    # --- Mapping CSV ---
    mapping_cols = [
        "target_taxonomy_id",
        "parent_taxonomy_id",
        "target_protein_uniprot_id",
        "target_protein_name",
        "target_protein_length",
        "target_organism",
        "support_pdb_id",
        "support_pdb_chain_description",
        "identity_pct",
        "alignment_score",
    ]

    if aggregated:
        df_map = pd.DataFrame(aggregated)[mapping_cols]
    else:
        df_map = pd.DataFrame(columns=mapping_cols)

    df_map.to_csv(mapping_path, index=False, encoding="utf-8-sig")
    log(f"输出映射文件: {mapping_path} ({len(df_map)} 行)")

    # --- All Alignments CSV ---
    aln_cols = [
        "target_accession",
        "target_protein_name",
        "target_length",
        "pdb_id",
        "chain_desc",
        "pdb_length",
        "identity_pct",
        "alignment_score",
    ]

    if all_alignments:
        df_aln = pd.DataFrame(all_alignments)[aln_cols]
    else:
        df_aln = pd.DataFrame(columns=aln_cols)

    df_aln.to_csv(alignments_path, index=False, encoding="utf-8-sig")
    log(f"输出比对详情: {alignments_path} ({len(df_aln)} 行)")

    # --- Log file ---
    with open(log_path, "w", encoding="utf-8") as f:
        f.write(f"Protein Homology Mapping Pipeline Log\n")
        f.write(f"{'='*50}\n")
        f.write(f"Target Taxonomy ID: {target_taxid}\n")
        f.write(f"Parent Taxonomy ID: {parent_taxid}\n")
        f.write(f"Target Organism: {target_name}\n")
        f.write(f"Identity Threshold: {thresholds}%\n")
        f.write(f"\n")
        f.write(f"Target Proteins (UniProt): {len(aggregated)} with matches\n")
        f.write(f"Total Alignments Passing Threshold: {len(all_alignments)}\n")
        f.write(f"\n")
        if aggregated:
            f.write(f"Top 10 Matches:\n")
            for row in aggregated[:10]:
                f.write(
                    f"  {row['target_protein_uniprot_id']}: "
                    f"{row['identity_pct']:.1f}% — {row['support_pdb_id']}\n"
                )

    log(f"输出日志文件: {log_path}")

    return mapping_path, alignments_path, log_path
