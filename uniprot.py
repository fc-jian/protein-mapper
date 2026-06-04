"""UniProt REST API client for fetching proteome sequences."""

from __future__ import annotations

import re
from typing import Any

import config
from utils import http_get, log


def fetch_uniprot_proteins(taxid: int) -> list[dict[str, Any]]:
    """Fetch all UniProt entries for a given taxonomy ID using paginated JSON search.

    Args:
        taxid: NCBI taxonomy ID.

    Returns:
        List of dicts, each containing:
            - accession: primaryAccession
            - entry_name: uniProtkbId
            - protein_name: recommended full name
            - sequence: amino acid sequence
            - length: sequence length
            - organism: scientific name
            - taxon_id: organism taxon ID
    """
    base_url = (
        f"https://rest.uniprot.org/uniprotkb/search"
        f"?query=taxonomy_id:{taxid}"
        f"&format=json"
        f"&size={config.UNIPROT_PAGE_SIZE}"
    )

    log(f"下载 UniProt 蛋白质组: taxonomy_id={taxid}")
    all_results: list[dict[str, Any]] = []
    url: str | None = base_url
    page = 0

    while url:
        page += 1
        log(f"  获取第 {page} 页...")
        resp = http_get(url)
        data = resp.json()

        results = data.get("results", [])
        all_results.extend(results)
        log(f"    本页 {len(results)} 条, 累计 {len(all_results)} 条")

        # Parse Link header for next page
        url = None
        link_header = resp.headers.get("Link", "")
        if link_header:
            for part in link_header.split(","):
                if 'rel="next"' in part:
                    url = part.split(";")[0].strip(" <>")
                    break

    proteins = _parse_uniprot_results(all_results)
    log(f"✓ UniProt 下载完成: 共 {len(proteins)} 个蛋白质")
    return proteins


def _parse_uniprot_results(results: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Parse raw UniProt JSON results into a clean protein dict list.

    Args:
        results: Raw list of UniProt JSON result entries.

    Returns:
        List of cleaned protein dicts.
    """
    proteins: list[dict[str, Any]] = []
    for entry in results:
        try:
            # Extract protein name safely
            protein_name = ""
            desc = entry.get("proteinDescription", {})
            recommended = desc.get("recommendedName", {})
            if recommended:
                protein_name = recommended.get("fullName", {}).get("value", "")
            if not protein_name:
                # Fallback to submitted names
                submitted = desc.get("submissionNames", [])
                if submitted:
                    protein_name = submitted[0].get("fullName", {}).get("value", "")

            seq_info = entry.get("sequence", {})
            organism_info = entry.get("organism", {})

            proteins.append({
                "accession": entry.get("primaryAccession", ""),
                "entry_name": entry.get("uniProtkbId", ""),
                "protein_name": protein_name,
                "sequence": seq_info.get("value", ""),
                "length": seq_info.get("length", 0),
                "organism": organism_info.get("scientificName", ""),
                "taxon_id": organism_info.get("taxonId", 0),
            })
        except Exception as exc:
            log(f"    解析条目失败: {exc}", "WARNING")
            continue

    return proteins
