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
from datetime import datetime

import config
from taxonomy import verify_lineage
from uniprot import fetch_uniprot_proteins
from pdb_search import search_pdb_structures, download_all_pdb_fasta
from alignment import run_alignment
from output import write_csv
from utils import log, probe_endpoints


def main() -> None:
    parser = argparse.ArgumentParser(
        description="蛋白质序列同源映射管道 — 目标物种蛋白组 × PDB 同源结构交叉比对",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--target", type=int, required=True, help="Target taxonomy ID")
    parser.add_argument("--parent", type=int, required=True, help="Parent taxonomy ID")
    parser.add_argument("--keywords", type=str, nargs="+", required=True, help="Space-separated PDB search keywords (AND logic)")
    parser.add_argument("--threshold", type=float, default=config.DEFAULT_THRESHOLD, help=f"Identity threshold %% (default: {config.DEFAULT_THRESHOLD})")
    parser.add_argument("--output", type=str, default="./", help="Output directory (default: ./)")
    args = parser.parse_args()

    target_taxid: int = args.target
    parent_taxid: int = args.parent
    keywords: list[str] = args.keywords
    threshold: float = args.threshold
    output_dir: str = args.output

    start_time = time.time()
    log("=" * 60)
    log("蛋白质序列同源映射管道 启动")
    log(f"  Target Taxonomy ID: {target_taxid}")
    log(f"  Parent Taxonomy ID: {parent_taxid}")
    log(f"  Keywords: {keywords}")
    log(f"  Identity Threshold: {threshold}%")
    log(f"  Output Directory: {output_dir}")
    log("=" * 60)

    # Probe all endpoint connectivity (cached for rest of pipeline)
    probe_endpoints()

    # 1. Verify taxonomy lineage
    log("\n[Step 1/5] 验证 NCBI Taxonomy 血缘关系")
    try:
        lineage_info = verify_lineage(target_taxid, parent_taxid)
        target_name = lineage_info["target_name"]
    except Exception as exc:
        log(f"  ✗ Lineage 验证失败: {exc}", "ERROR")
        sys.exit(1)

    # 2. Download target proteome from UniProt
    log("\n[Step 2/5] 下载目标物种蛋白质组 (UniProt)")
    try:
        target_proteins = fetch_uniprot_proteins(target_taxid)
    except Exception as exc:
        log(f"  ✗ UniProt 下载失败: {exc}", "ERROR")
        sys.exit(1)

    if not target_proteins:
        log("  ⚠ UniProt 未返回任何蛋白质, 退出", "WARNING")
        write_csv([], [], target_taxid, parent_taxid, threshold, output_dir, target_name)
        sys.exit(0)

    # 3. Search PDB structures
    log("\n[Step 3/5] 搜索 RCSB PDB 结构")
    try:
        pdb_ids = search_pdb_structures(parent_taxid, keywords)
    except Exception as exc:
        log(f"  ✗ PDB 搜索失败: {exc}", "ERROR")
        sys.exit(1)

    if not pdb_ids:
        log("  ⚠ 未找到匹配的 PDB 结构, 退出", "WARNING")
        write_csv([], [], target_taxid, parent_taxid, threshold, output_dir, target_name)
        sys.exit(0)

    # 4. Download PDB FASTA
    log("\n[Step 4/5] 下载 PDB FASTA 序列")
    try:
        pdb_chains = download_all_pdb_fasta(pdb_ids)
    except Exception as exc:
        log(f"  ✗ PDB FASTA 下载失败: {exc}", "ERROR")
        sys.exit(1)

    if not pdb_chains:
        log("  ⚠ 未获取到 PDB chain 序列, 退出", "WARNING")
        write_csv([], [], target_taxid, parent_taxid, threshold, output_dir, target_name)
        sys.exit(0)

    # 5. Sequence alignment
    log("\n[Step 5/5] 序列比对")
    aggregated, all_alignments = run_alignment(target_proteins, pdb_chains, threshold)

    # 7. Write output
    log("\n[输出] 写入 CSV 文件")
    mapping_path, aln_path, log_path = write_csv(
        aggregated, all_alignments, target_taxid, parent_taxid, threshold, output_dir, target_name
    )

    elapsed = time.time() - start_time
    log("=" * 60)
    log(f"✓ 管道执行完成! 耗时 {elapsed:.1f}s")
    log(f"  映射文件: {mapping_path}")
    log(f"  比对详情: {aln_path}")
    log(f"  运行日志: {log_path}")
    log(f"  统计: {len(target_proteins)} 目标蛋白, {len(pdb_ids)} PDB, {len(aggregated)} 蛋白有匹配")
    log("=" * 60)


if __name__ == "__main__":
    main()
