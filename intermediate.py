"""Writers for intermediate sequence and metadata files."""

from __future__ import annotations

import csv
import json
import os
from pathlib import Path
from typing import Any

import config
from utils import log


def intermediate_root(output_dir: str, target_taxid: int, parent_taxid: int) -> Path:
    """Return the default intermediate directory for one pipeline run."""
    base = f"{target_taxid}_{parent_taxid}{config.INTERMEDIATE_DIR_SUFFIX}"
    return Path(output_dir) / base


def write_fasta_record(handle: Any, record_id: str, sequence: str, description: str = "") -> None:
    """Write one FASTA record with wrapped sequence lines."""
    header = record_id if not description else f"{record_id} {description}"
    handle.write(f">{header}\n")
    for start in range(0, len(sequence), 80):
        handle.write(f"{sequence[start:start + 80]}\n")


def save_uniprot_intermediates(
    proteins: list[dict[str, Any]],
    output_dir: str,
    target_taxid: int,
    parent_taxid: int,
) -> tuple[str, str]:
    """Save normalized UniProt proteins as FASTA and CSV metadata."""
    if not config.SAVE_INTERMEDIATES:
        return "", ""

    out_dir = intermediate_root(output_dir, target_taxid, parent_taxid) / "uniprot"
    out_dir.mkdir(parents=True, exist_ok=True)
    fasta_path = out_dir / f"{target_taxid}_uniprot.fasta"
    metadata_path = out_dir / f"{target_taxid}_uniprot_proteins.csv"

    with fasta_path.open("w", encoding="utf-8") as handle:
        for protein in proteins:
            description = " ".join(
                part
                for part in [
                    protein.get("entry_name", ""),
                    protein.get("protein_name", ""),
                    f"OS={protein.get('organism', '')}" if protein.get("organism") else "",
                    f"OX={protein.get('taxon_id', '')}" if protein.get("taxon_id") else "",
                ]
                if part
            )
            write_fasta_record(handle, protein.get("accession", ""), protein.get("sequence", ""), description)

    fieldnames = ["accession", "entry_name", "protein_name", "length", "organism", "taxon_id"]
    with metadata_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for protein in proteins:
            writer.writerow({key: protein.get(key, "") for key in fieldnames})

    log(f"Saved UniProt intermediate FASTA: {fasta_path}")
    log(f"Saved UniProt intermediate metadata: {metadata_path}")
    return str(fasta_path), str(metadata_path)


def save_pdb_search_intermediates(
    pdb_ids: list[str],
    keywords: list[str],
    output_dir: str,
    target_taxid: int,
    parent_taxid: int,
) -> tuple[str, str]:
    """Save PDB search IDs and query metadata."""
    if not config.SAVE_INTERMEDIATES:
        return "", ""

    out_dir = intermediate_root(output_dir, target_taxid, parent_taxid) / "pdb"
    out_dir.mkdir(parents=True, exist_ok=True)
    ids_path = out_dir / f"{target_taxid}_{parent_taxid}_pdb_ids.txt"
    query_path = out_dir / f"{target_taxid}_{parent_taxid}_pdb_query.json"

    ids_path.write_text("\n".join(pdb_ids) + ("\n" if pdb_ids else ""), encoding="utf-8")
    query_path.write_text(
        json.dumps(
            {
                "target_taxonomy_id": target_taxid,
                "parent_taxonomy_id": parent_taxid,
                "keywords": keywords,
                "pdb_count": len(pdb_ids),
            },
            indent=2,
            ensure_ascii=False,
        ),
        encoding="utf-8",
    )

    log(f"Saved PDB ID list: {ids_path}")
    log(f"Saved PDB query metadata: {query_path}")
    return str(ids_path), str(query_path)


def pdb_fasta_dir(output_dir: str, target_taxid: int, parent_taxid: int) -> str:
    """Return the default directory for raw downloaded PDB FASTA files."""
    return str(intermediate_root(output_dir, target_taxid, parent_taxid) / "pdb" / "fasta")


def save_pdb_chain_intermediates(
    pdb_chains: list[dict[str, Any]],
    output_dir: str,
    target_taxid: int,
    parent_taxid: int,
) -> str:
    """Save parsed PDB chain metadata as CSV."""
    if not config.SAVE_INTERMEDIATES:
        return ""

    out_dir = intermediate_root(output_dir, target_taxid, parent_taxid) / "pdb"
    out_dir.mkdir(parents=True, exist_ok=True)
    chains_path = out_dir / f"{target_taxid}_{parent_taxid}_pdb_chains.csv"
    fieldnames = ["pdb_id", "chain_desc", "length", "sequence"]

    with chains_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for chain in pdb_chains:
            writer.writerow({key: chain.get(key, "") for key in fieldnames})

    log(f"Saved PDB chain metadata: {chains_path}")
    return str(chains_path)


def save_raw_pdb_fasta(fasta_dir: str | None, pdb_id: str, fasta_text: str) -> str:
    """Save one raw PDB FASTA response if an output directory is configured."""
    if not fasta_dir:
        return ""

    os.makedirs(fasta_dir, exist_ok=True)
    path = Path(fasta_dir) / f"{pdb_id}.fasta"
    path.write_text(fasta_text, encoding="utf-8")
    return str(path)
