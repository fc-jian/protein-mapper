"""MMseqs2 clustering for mapped UniProt proteins."""

from __future__ import annotations

import csv
import shutil
import subprocess
from pathlib import Path
from typing import Any

import pandas as pd

import config
from intermediate import intermediate_root, write_fasta_record
from output import MAPPING_COLUMNS
from utils import log


CLUSTER_COLUMNS = [
    "cluster_id",
    "cluster_representative",
    "cluster_members",
    "cluster_size",
    *MAPPING_COLUMNS,
]


def run_mmseqs_clustering_report(
    aggregated: list[dict[str, Any]],
    target_proteins: list[dict[str, Any]],
    target_taxid: int,
    parent_taxid: int,
    output_dir: str,
    cluster_threshold: float = config.MMSEQS_CLUSTER_THRESHOLD,
) -> str:
    """Cluster mapped UniProt sequences and write a cluster-level report CSV."""
    base = f"{target_taxid}_{parent_taxid}"
    cluster_csv_path = Path(output_dir) / f"{base}_cluster_mapping.csv"

    if not aggregated:
        _write_cluster_csv([], cluster_csv_path)
        log(f"No mapped proteins to cluster; wrote empty cluster report: {cluster_csv_path}")
        return str(cluster_csv_path)

    protein_by_accession = {protein["accession"]: protein for protein in target_proteins}
    mapped_accessions = [row["target_protein_uniprot_id"] for row in aggregated]
    mapped_proteins = [
        protein_by_accession[accession]
        for accession in mapped_accessions
        if accession in protein_by_accession and protein_by_accession[accession].get("sequence")
    ]

    if not mapped_proteins:
        _write_cluster_csv([], cluster_csv_path)
        log(f"No mapped UniProt sequences available for clustering; wrote empty cluster report: {cluster_csv_path}")
        return str(cluster_csv_path)

    mmseqs_dir = intermediate_root(output_dir, target_taxid, parent_taxid) / "mmseqs"
    mmseqs_dir.mkdir(parents=True, exist_ok=True)
    input_fasta = mmseqs_dir / f"{base}_mapped_uniprot.fasta"
    output_prefix = mmseqs_dir / f"{base}_mmseqs_cluster"
    tmp_dir = mmseqs_dir / "tmp"
    cluster_tsv = Path(f"{output_prefix}_cluster.tsv")

    _write_mapped_fasta(mapped_proteins, input_fasta)
    _clear_previous_mmseqs_outputs(mmseqs_dir, output_prefix.name, tmp_dir)

    log(
        f"Clustering {len(mapped_proteins)} mapped UniProt sequences with MMseqs2 "
        f"(min sequence identity={cluster_threshold})"
    )
    command = [
        config.MMSEQS_BINARY,
        "easy-cluster",
        str(input_fasta),
        str(output_prefix),
        str(tmp_dir),
        "--min-seq-id",
        str(cluster_threshold),
        "--threads",
        str(config.MMSEQS_THREADS),
    ]

    try:
        completed = subprocess.run(command, check=False, capture_output=True, text=True)
    except FileNotFoundError as exc:
        raise RuntimeError(
            f"MMseqs2 executable '{config.MMSEQS_BINARY}' was not found. "
            "Install MMseqs2 or adjust config.MMSEQS_BINARY."
        ) from exc

    if completed.returncode != 0:
        stderr = completed.stderr.strip() or completed.stdout.strip()
        raise RuntimeError(f"MMseqs2 clustering failed with exit code {completed.returncode}: {stderr[-2000:]}")

    clusters = _read_mmseqs_clusters(cluster_tsv)
    cluster_rows = _build_cluster_rows(clusters, aggregated)
    _write_cluster_csv(cluster_rows, cluster_csv_path)
    log(f"Wrote cluster-level mapping report: {cluster_csv_path} ({len(cluster_rows)} rows)")
    return str(cluster_csv_path)


def _write_mapped_fasta(proteins: list[dict[str, Any]], path: Path) -> None:
    with path.open("w", encoding="utf-8") as handle:
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
            write_fasta_record(handle, protein["accession"], protein["sequence"], description)
    log(f"Saved mapped UniProt FASTA for MMseqs2: {path}")


def _clear_previous_mmseqs_outputs(mmseqs_dir: Path, output_prefix_name: str, tmp_dir: Path) -> None:
    for path in mmseqs_dir.glob(f"{output_prefix_name}*"):
        if path.is_dir():
            shutil.rmtree(path)
        else:
            path.unlink()
    if tmp_dir.exists():
        shutil.rmtree(tmp_dir)


def _read_mmseqs_clusters(cluster_tsv: Path) -> dict[str, list[str]]:
    clusters: dict[str, list[str]] = {}
    with cluster_tsv.open("r", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if len(row) < 2:
                continue
            representative, member = row[0], row[1]
            clusters.setdefault(representative, [])
            if member not in clusters[representative]:
                clusters[representative].append(member)
    return clusters


def _build_cluster_rows(
    clusters: dict[str, list[str]],
    aggregated: list[dict[str, Any]],
) -> list[dict[str, str]]:
    rows_by_accession = {row["target_protein_uniprot_id"]: row for row in aggregated}
    pending_rows: list[dict[str, str]] = []

    for representative, members in clusters.items():
        known_members = [member for member in members if member in rows_by_accession]
        if not known_members:
            continue

        member_rows = [rows_by_accession[member] for member in known_members]
        cluster_row: dict[str, str] = {
            "cluster_id": "",
            "cluster_representative": representative,
            "cluster_members": ";".join(known_members),
            "cluster_size": str(len(known_members)),
        }

        for column in MAPPING_COLUMNS:
            cluster_row[column] = _dedup_join(row.get(column, "") for row in member_rows)

        pending_rows.append(cluster_row)

    pending_rows.sort(key=lambda row: int(row["cluster_size"]), reverse=True)
    for index, row in enumerate(pending_rows, start=1):
        row["cluster_id"] = f"cluster_{index}"
    return pending_rows


def _dedup_join(values: Any) -> str:
    seen: set[str] = set()
    joined: list[str] = []
    for value in values:
        for part in str(value).split(";"):
            cleaned = part.strip()
            if cleaned and cleaned not in seen:
                seen.add(cleaned)
                joined.append(cleaned)
    return ";".join(joined)


def _write_cluster_csv(rows: list[dict[str, str]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if rows:
        df = pd.DataFrame(rows)[CLUSTER_COLUMNS]
    else:
        df = pd.DataFrame(columns=CLUSTER_COLUMNS)
    df.to_csv(path, index=False, encoding="utf-8-sig")
