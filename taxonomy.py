"""NCBI Taxonomy lineage verification using the Datasets v2 API."""

from __future__ import annotations

from typing import Any

from utils import http_get, log


def verify_lineage(target_taxid: int, parent_taxid: int) -> dict[str, Any]:
    """Verify that parent_taxid is an ancestor of target_taxid via NCBI lineage.

    Fetches the taxonomy node for target_taxid and checks whether parent_taxid
    appears in its lineage array.

    Args:
        target_taxid: The target taxonomy ID (e.g. 10253 for Vaccinia virus Tian Tan).
        parent_taxid: The parent taxonomy ID to verify (e.g. 10242 for Chordopoxvirinae).

    Returns:
        Dict with 'target_name' and 'lineage' keys.

    Raises:
        ValueError: If parent_taxid is not in the target's lineage.
        requests.RequestException: On API failure.
    """
    url = f"https://api.ncbi.nlm.nih.gov/datasets/v2/taxonomy/taxon/{target_taxid}"
    log(f"Verifying NCBI lineage: target={target_taxid}, parent={parent_taxid}")

    resp = http_get(url)
    data = resp.json()

    nodes = data.get("taxonomy_nodes", [])
    if not nodes:
        raise ValueError(f"NCBI returned no taxonomy node for tax_id={target_taxid}")

    taxonomy = nodes[0].get("taxonomy", {})
    organism_name = taxonomy.get("organism_name", str(target_taxid))
    lineage: list[int] = taxonomy.get("lineage", [])

    log(f"  Target organism: {organism_name}")
    log(f"  Lineage preview (first 5 + last 3): {lineage[:5]}...{lineage[-3:] if len(lineage) > 8 else ''}")

    if parent_taxid not in lineage:
        raise ValueError(
            f"Parent {parent_taxid} is not an ancestor of target {target_taxid} ({organism_name}). "
            f"{parent_taxid} was not found in the lineage."
        )

    log(f"  Lineage verification passed: {parent_taxid} is an ancestor of {target_taxid}")
    return {"target_name": organism_name, "lineage": lineage}
