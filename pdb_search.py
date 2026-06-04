"""RCSB PDB Search API v2 client — structure search + FASTA download."""

from __future__ import annotations

import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from io import StringIO
from typing import Any

from Bio import SeqIO

import config
from utils import http_get, http_post, log


def search_pdb_structures(
    parent_taxid: int,
    keywords: list[str],
) -> list[str]:
    """Search RCSB PDB for structures matching taxonomy lineage + keywords.

    Constructs a compound query with:
      - taxonomy_lineage.id exact_match → parent_taxid
      - (struct.title OR rcsb_polymer_entity.pdbx_description) contains_words → keywords (any-of / OR logic)

    Args:
        parent_taxid: Parent taxonomy ID for lineage filtering.
        keywords: List of keywords (AND logic via space-joining).

    Returns:
        List of PDB IDs (e.g. ['2I9L', '5EOR']).
    """
    keyword_str = " ".join(keywords)
    log(f"搜索 RCSB PDB: parent_taxonomy={parent_taxid}, keywords={keyword_str}")

    # Build OR nodes: one for struct.title, one for polymer entity description
    keyword_nodes = []
    for attr in ["struct.title", "rcsb_polymer_entity.pdbx_description"]:
        keyword_nodes.append({
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": attr,
                "operator": "contains_words",
                "value": keyword_str,
            },
        })

    query = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entity_source_organism.taxonomy_lineage.id",
                        "operator": "exact_match",
                        "value": str(parent_taxid),
                    },
                },
                {
                    "type": "group",
                    "logical_operator": "or",
                    "nodes": keyword_nodes,
                },
            ],
        },
        "return_type": "entry",
        "request_options": {
            "paginate": {"start": 0, "rows": config.PDB_PAGE_SIZE},
        },
    }

    all_ids: list[str] = []
    start = 0

    while True:
        query["request_options"]["paginate"]["start"] = start  # type: ignore[index]
        resp = http_post("https://search.rcsb.org/rcsbsearch/v2/query", query)
        data = resp.json()

        total = data.get("total_count", 0)
        result_set = data.get("result_set", [])
        batch_ids = [item["identifier"] for item in result_set]
        all_ids.extend(batch_ids)

        log(f"  获取 {len(batch_ids)} 个 PDB (start={start}, total={total})")

        start += len(batch_ids)
        if start >= total or len(batch_ids) == 0:
            break

    log(f"✓ PDB 搜索完成: 共 {len(all_ids)} 个结构")
    return all_ids


def fetch_pdb_fasta(pdb_id: str) -> tuple[str, str | None]:
    """Download FASTA for a single PDB entry.

    Args:
        pdb_id: PDB identifier (e.g. '2I9L').

    Returns:
        Tuple of (pdb_id, fasta_text). fasta_text is None on failure.
    """
    url = f"https://www.rcsb.org/fasta/entry/{pdb_id}/download"
    try:
        resp = http_get(url, timeout=config.FASTA_TIMEOUT, max_retries=2)
        return pdb_id, resp.text
    except Exception as exc:
        log(f"  PDB {pdb_id} FASTA 下载失败: {exc}", "WARNING")
        return pdb_id, None


def download_all_pdb_fasta(
    pdb_ids: list[str],
) -> list[dict[str, Any]]:
    """Download and parse FASTA for multiple PDB entries concurrently.

    Args:
        pdb_ids: List of PDB identifiers.

    Returns:
        List of chain dicts with keys: pdb_id, chain_desc, sequence, length.
    """
    if not pdb_ids:
        log("  没有 PDB ID 需要下载", "WARNING")
        return []

    log(f"下载 {len(pdb_ids)} 个 PDB 的 FASTA 序列 (并发数={config.PDB_CONCURRENT_WORKERS})...")

    all_chains: list[dict[str, Any]] = []
    success = 0
    failed = 0

    with ThreadPoolExecutor(max_workers=config.PDB_CONCURRENT_WORKERS) as executor:
        futures = {
            executor.submit(fetch_pdb_fasta, pid): pid
            for pid in pdb_ids
        }
        for future in as_completed(futures):
            pdb_id, fasta_text = future.result()
            if fasta_text is None:
                failed += 1
                continue
            success += 1
            chains = _parse_pdb_fasta(pdb_id, fasta_text)
            all_chains.extend(chains)

    log(f"✓ PDB FASTA 下载完成: 成功 {success}, 失败 {failed}, 共 {len(all_chains)} 条 chain")
    return all_chains


def _parse_pdb_fasta(pdb_id: str, fasta_text: str) -> list[dict[str, Any]]:
    """Parse a multi-sequence PDB FASTA into individual chain records.

    FASTA header format:
        >{PDB_ID}_{entity_id}|Chain(s) {chain_ids}|{protein_name}|{organism} ({tax_id})

    Args:
        pdb_id: The PDB identifier.
        fasta_text: Raw FASTA text from RCSB.

    Returns:
        List of chain dicts with: pdb_id, chain_desc, sequence, length.
    """
    chains: list[dict[str, Any]] = []
    try:
        for record in SeqIO.parse(StringIO(fasta_text), "fasta"):
            seq = str(record.seq)
            if seq:
                chains.append({
                    "pdb_id": pdb_id,
                    "chain_desc": record.description,
                    "sequence": seq,
                    "length": len(seq),
                })
    except Exception as exc:
        log(f"  解析 PDB {pdb_id} FASTA 失败: {exc}", "WARNING")

    return chains
