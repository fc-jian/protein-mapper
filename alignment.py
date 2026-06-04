"""Sequence alignment engine using BioPython PairwiseAligner with k-mer pre-screening."""

from __future__ import annotations

from typing import Any

from Bio.Align import PairwiseAligner

import config
from utils import log

# Reusable aligner instance configured for local alignment
_aligner = PairwiseAligner()
_aligner.mode = "local"
_aligner.open_gap_score = -10.0
_aligner.extend_gap_score = -0.5


def kmer_jaccard(seq1: str, seq2: str, k: int = config.KMER_SIZE) -> float:
    """Compute k-mer Jaccard similarity for fast pre-screening.

    Args:
        seq1: First amino acid sequence.
        seq2: Second amino acid sequence.
        k: k-mer size (default 3 for tripeptides).

    Returns:
        Jaccard index in [0, 1].
    """
    if len(seq1) < k or len(seq2) < k:
        return 0.0

    kmers1 = {seq1[i : i + k] for i in range(len(seq1) - k + 1)}
    kmers2 = {seq2[i : i + k] for i in range(len(seq2) - k + 1)}

    intersection = len(kmers1 & kmers2)
    union = len(kmers1 | kmers2)
    return intersection / union if union > 0 else 0.0


def align_and_identity(seq1: str, seq2: str) -> tuple[float, float]:
    """Perform local pairwise alignment and compute sequence identity.

    Uses Bio.Align.PairwiseAligner with local mode and configured gap scores.

    Args:
        seq1: Query amino acid sequence.
        seq2: Subject amino acid sequence.

    Returns:
        Tuple of (identity_percent, alignment_score). Identity is 0.0 if no
        alignment is found.
    """
    alignments = _aligner.align(seq1, seq2)
    if not alignments:
        return 0.0, 0.0

    # Take the best alignment
    best = alignments[0]
    aligned1 = str(best[0])
    aligned2 = str(best[1])
    score = best.score
    aligned_len = len(aligned1)
    if aligned_len == 0:
        return 0.0, 0.0

    matches = sum(1 for a, b in zip(aligned1, aligned2) if a == b)
    identity = (matches / aligned_len) * 100.0
    return identity, score


def run_alignment(
    target_proteins: list[dict[str, Any]],
    pdb_chains: list[dict[str, Any]],
    threshold: float = config.DEFAULT_THRESHOLD,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Run all-vs-all alignment between target proteins and PDB chains.

    Applies length-ratio pre-filter, optional k-mer pre-screening for large
    datasets, and pairwise alignment. Only results with identity >= threshold
    are kept.

    Args:
        target_proteins: List of target protein dicts from UniProt.
        pdb_chains: List of PDB chain dicts from RCSB.
        threshold: Minimum identity percentage to keep a match.

    Returns:
        Tuple of:
            - aggregated_results: one row per target protein (for summary CSV).
            - all_alignments: every qualifying alignment (for detail CSV).
    """
    n_proteins = len(target_proteins)
    n_chains = len(pdb_chains)
    total_alignments = n_proteins * n_chains

    log(f"Starting sequence alignment: {n_proteins} target proteins x {n_chains} PDB chains = {total_alignments} comparisons")
    log(f"  Identity threshold: {threshold}%")

    use_kmer = total_alignments > config.KMER_PRE_FILTER_TRIGGER
    if use_kmer:
        log(f"  Using k-mer pre-screening (k={config.KMER_SIZE}), Jaccard >= {config.KMER_JACCARD_THRESHOLD}")

    all_alignments: list[dict[str, Any]] = []
    matched_count = 0
    skipped_length = 0
    skipped_kmer = 0

    for i, protein in enumerate(target_proteins):
        seq1 = protein["sequence"]
        len1 = protein["length"]
        if i % 50 == 0 and i > 0:
            log(f"  Alignment progress: {i}/{n_proteins} proteins, {matched_count} matches")

        for chain in pdb_chains:
            seq2 = chain["sequence"]
            len2 = chain["length"]

            # Length ratio pre-filter
            if len1 > 0 and len2 > 0:
                ratio = max(len1, len2) / min(len1, len2)
                if ratio > config.LENGTH_RATIO_MAX:
                    skipped_length += 1
                    continue

            # k-mer pre-screening
            if use_kmer:
                jaccard = kmer_jaccard(seq1, seq2)
                if jaccard < config.KMER_JACCARD_THRESHOLD:
                    skipped_kmer += 1
                    continue

            # Pairwise alignment
            identity, score = align_and_identity(seq1, seq2)
            if identity >= threshold:
                alignment_record = {
                    "target_accession": protein["accession"],
                    "target_entry_name": protein["entry_name"],
                    "target_protein_name": protein["protein_name"],
                    "target_sequence": seq1,
                    "target_length": len1,
                    "pdb_id": chain["pdb_id"],
                    "chain_desc": chain["chain_desc"],
                    "pdb_sequence": seq2,
                    "pdb_length": len2,
                    "identity_pct": round(identity, 2),
                    "alignment_score": round(score, 2),
                }
                all_alignments.append(alignment_record)
                matched_count += 1

    log(f"  Alignment stats: length-filtered {skipped_length}, k-mer-filtered {skipped_kmer}, matches {matched_count}")

    # Aggregate by target protein: best chain per PDB
    aggregated = _aggregate_results(all_alignments, target_proteins)
    log(f"Alignment completed: {len(aggregated)} target proteins have PDB matches ({len(all_alignments)} alignment records)")
    return aggregated, all_alignments


def _aggregate_results(
    all_alignments: list[dict[str, Any]],
    target_proteins: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Aggregate alignment results: one row per target protein.

    For each target protein, groups matches by PDB ID and keeps the
    highest-identity chain per PDB.

    Args:
        all_alignments: All qualifying alignment records.
        target_proteins: Complete list of target proteins (to include those
            with no matches).

    Returns:
        Aggregated list, one dict per target protein.
    """
    # Group alignments by target accession
    by_accession: dict[str, list[dict[str, Any]]] = {}
    for aln in all_alignments:
        acc = aln["target_accession"]
        by_accession.setdefault(acc, []).append(aln)

    # Build target protein lookup
    protein_lookup: dict[str, dict[str, Any]] = {
        p["accession"]: p for p in target_proteins
    }

    aggregated: list[dict[str, Any]] = []
    for acc, alignments in by_accession.items():
        # Best chain per PDB
        best_per_pdb: dict[str, tuple[str, float, float]] = {}
        for aln in alignments:
            pid = aln["pdb_id"]
            cdesc = aln["chain_desc"]
            ident = aln["identity_pct"]
            score = aln["alignment_score"]
            if pid not in best_per_pdb or ident > best_per_pdb[pid][1]:
                best_per_pdb[pid] = (cdesc, ident, score)

        pdb_ids = sorted(best_per_pdb.keys())
        protein = protein_lookup.get(acc, {})

        aggregated.append({
            "target_taxonomy_id": protein.get("taxon_id", ""),
            "parent_taxonomy_id": "",  # filled by caller
            "target_protein_uniprot_id": acc,
            "target_protein_name": protein.get("protein_name", ""),
            "target_protein_length": protein.get("length", 0),
            "target_organism": protein.get("organism", ""),
            "support_pdb_id": ";".join(pdb_ids),
            "support_pdb_chain_description": ";".join(best_per_pdb[pid][0] for pid in pdb_ids),
            "identity_pct": max(v[1] for v in best_per_pdb.values()),
            "alignment_score": max(v[2] for v in best_per_pdb.values()),
        })

    # Sort by identity descending
    aggregated.sort(key=lambda x: x["identity_pct"], reverse=True)
    return aggregated
