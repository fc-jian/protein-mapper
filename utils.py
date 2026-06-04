"""HTTP request utilities with endpoint probes and retry logic."""

from __future__ import annotations

import time
from datetime import datetime
from typing import Any

import requests

import config


# ── endpoint connectivity probe ─────────────────────────────────────────────

_probe_cache: dict[str, bool] = {}


def probe_endpoints() -> None:
    """Probe all known endpoints once at startup."""
    endpoints = {
        "ncbi": "https://api.ncbi.nlm.nih.gov/datasets/v2/taxonomy/taxon/1",
        "uniprot": "https://rest.uniprot.org/uniprotkb/search?query=taxonomy_id:1&format=json&size=1",
        "rcsb_search": "https://search.rcsb.org/rcsbsearch/v2/query",
        "rcsb_fasta": "https://www.rcsb.org/fasta/entry/1AKE/download",
    }

    for name, url in endpoints.items():
        if _probe_cache.get(name):
            continue

        log(f"探测端点连通性 [{name}]: {url[:70]}...")
        if name == "rcsb_search":
            _probe_one_post(
                url,
                {
                    "query": {
                        "type": "terminal",
                        "service": "text",
                        "parameters": {
                            "attribute": "struct.title",
                            "operator": "contains_phrase",
                            "value": "probe",
                        },
                    },
                    "return_type": "entry",
                    "request_options": {"paginate": {"start": 0, "rows": 1}},
                },
            )
        else:
            _probe_one_get(url)

        _probe_cache[name] = True
        log(f"  端点 [{name}]: 可用")

    status = ", ".join(f"{name}=可用" for name in _probe_cache)
    log(f"端点探测完成: {status}")


def _probe_one_get(url: str) -> None:
    """Probe one GET endpoint and raise RuntimeError on failure."""
    try:
        resp = requests.get(url, timeout=config.PROBE_TIMEOUT)
        resp.raise_for_status()
    except requests.RequestException as exc:
        log(f"  端点不可用: {url}", "ERROR")
        raise RuntimeError(f"无法连接 {url}") from exc


def _probe_one_post(url: str, json_data: dict[str, Any]) -> None:
    """Probe one POST endpoint and raise RuntimeError on failure."""
    try:
        resp = requests.post(url, json=json_data, timeout=config.PROBE_TIMEOUT)
        resp.raise_for_status()
    except requests.RequestException as exc:
        log(f"  端点不可用: {url}", "ERROR")
        raise RuntimeError(f"无法连接 {url}") from exc


# ── logging ─────────────────────────────────────────────────────────────────


def log(message: str, level: str = "INFO") -> None:
    """Print a timestamped log message."""
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{ts}] [{level}] {message}", flush=True)


# ── core request helpers ────────────────────────────────────────────────────


def _try_request(
    method: str,
    url: str,
    *,
    timeout: int = config.TIMEOUT,
    max_retries: int = config.MAX_RETRIES,
    json_data: dict[str, Any] | None = None,
) -> requests.Response:
    """Issue an HTTP request with retries."""
    last_exc: Exception | None = None

    for attempt in range(max_retries + 1):
        try:
            if method == "GET":
                resp = requests.get(url, timeout=timeout)
            else:
                resp = requests.post(url, json=json_data, timeout=timeout)
            resp.raise_for_status()
            return resp
        except requests.RequestException as exc:
            last_exc = exc
            if attempt < max_retries:
                wait = config.RETRY_BACKOFF_BASE ** (attempt + 1)
                log(
                    f"{method} {url[:80]}... 失败 "
                    f"(尝试 {attempt + 1}/{max_retries}): {exc!r:.100}。"
                    f"{wait:.0f}s 后重试...",
                    "WARNING",
                )
                time.sleep(wait)

    raise last_exc  # type: ignore[misc]


def http_get(
    url: str,
    timeout: int = config.TIMEOUT,
    max_retries: int = config.MAX_RETRIES,
) -> requests.Response:
    """HTTP GET with retry handling."""
    return _try_request("GET", url, timeout=timeout, max_retries=max_retries)


def http_post(
    url: str,
    json_data: dict[str, Any],
    timeout: int = config.TIMEOUT,
    max_retries: int = config.MAX_RETRIES,
) -> requests.Response:
    """HTTP POST with retry handling."""
    return _try_request(
        "POST",
        url,
        timeout=timeout,
        max_retries=max_retries,
        json_data=json_data,
    )
