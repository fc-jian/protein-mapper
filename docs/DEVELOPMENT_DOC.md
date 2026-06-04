# Protein Sequence Homology Mapping Pipeline — 开发文档

## 1. 概述

本程序自动化完成「目标物种蛋白质组 × PDB同源结构序列」的交叉比对任务：

```
输入: target_taxonomy_id, parent_taxonomy_id, add_keywords[]
  │
  ├─ Step 1: NCBI Taxonomy 验证 parent ⊃ target (lineage check)
  ├─ Step 2: UniProt → 下载 target 全蛋白质组的完整序列
  ├─ Step 3: RCSB PDB → 搜索 parent 下含所有 add_keywords 的 PDB 结构
  ├─ Step 4: RCSB PDB → 逐 PDB 下载 FASTA 序列 (含 chain 描述)
  └─ Step 5: 序列比对 → 输出 CSV (identity ≥ threshold)
```

示例输入:
- `target_taxonomy_id = 10253` (Vaccinia virus Tian Tan)
- `parent_taxonomy_id = 10242` (Chordopoxvirinae / 痘病毒亚科)
- `add_keywords = ["antibody"]`

逻辑意义: 在 parent 类群的所有 PDB 抗体结构中寻找与 target 同源蛋白有序列相似性的 chain，这些 chain 可能是 target 蛋白的抗原表位所在结构域。

---

## 2. 环境准备

### 2.1 网络访问

程序需要能访问 NCBI、UniProt 与 RCSB PDB 的公开 HTTP/HTTPS API。

### 2.2 Python 依赖

```
biopython>=1.83    # SeqIO (FASTA 解析), PairwiseAligner (序列比对)
requests>=2.31      # HTTP 请求
pandas>=2.0         # CSV 输出
```

安装:
```bash
pip install biopython requests pandas
```

---

## 3. API 详细规格

### 3.1 NCBI Datasets v2 — Taxonomy Lineage 验证

**端点**: `GET https://api.ncbi.nlm.nih.gov/datasets/v2/taxonomy/taxon/{tax_id}`

**响应结构** (JSON):
```json
{
  "taxonomy_nodes": [{
    "taxonomy": {
      "tax_id": 10253,
      "organism_name": "Vaccinia virus Tian Tan",
      "lineage": [1, 10239, ..., 10242, 10245]  // 从根到自身的祖先链
    }
  }]
}
```

**验证逻辑**:
```python
if parent_taxid not in lineage:
    raise ValueError(f"Parent {parent_taxid} is not an ancestor of target {target_taxid}")
```

**注意**: `lineage` 数组包含自身 ID，因此 parent 在 lineage 中即表示 target 是 parent 的后代。

---

### 3.2 UniProt REST API — 获取目标物种蛋白质组

**端点 (search)**: `GET https://rest.uniprot.org/uniprotkb/search`

| 参数 | 值 | 说明 |
|------|-----|------|
| `query` | `taxonomy_id:{TAXID}` | 按 taxonomy 过滤 |
| `format` | `json` | JSON 格式(含序列) |
| `size` | `500` | 每页最多 500 条 |

**分页**: 响应头 `Link` 包含下一页 URL (含 `cursor` 参数):
```
Link: <https://rest.uniprot.org/uniprotkb/search?...&cursor=xxx>; rel="next"
```

**端点 (FASTA stream)**:
`GET https://rest.uniprot.org/uniprotkb/stream?query=taxonomy_id:{TAXID}&format=fasta`

- 返回多序列 FASTA，含完整 header line
- FASTA header 格式: `>sp|ACCESSION|ENTRY_NAME Description OS=Organism OX=TAXID GN=GeneName PE=N SV=N`
- 流式返回，适合批量下载

**推荐策略**: 使用 **search (JSON)** 逐页获取，JSON 数据结构化程度高，同时包含 `primaryAccession`、`uniProtkbId`、`proteinDescription`、`sequence.value` 等字段，便于后续输出 CSV 中的 protein name。

**JSON 结果关键字段** (从 `results[]` 中提取):
| 字段路径 | 说明 |
|----------|------|
| `primaryAccession` | UniProt accession (如 Q9JFA8) |
| `uniProtkbId` | Entry name (如 MCEL_VACCT) |
| `proteinDescription.recommendedName.fullName.value` | 蛋白名 |
| `sequence.value` | 氨基酸序列 (单字母) |
| `sequence.length` | 序列长度 |
| `organism.scientificName` | 物种名 |
| `organism.taxonId` | Taxonomy ID |

**分页实现伪代码**:
```python
def fetch_all_uniprot(taxid):
    url = f"https://rest.uniprot.org/uniprotkb/search?query=taxonomy_id:{taxid}&format=json&size=500"
    all_results = []
    while url:
        resp = requests.get(url, timeout=60)
        resp.raise_for_status()
        data = resp.json()
        all_results.extend(data['results'])
        # Parse Link header for next page
        url = None
        if 'Link' in resp.headers:
            for part in resp.headers['Link'].split(','):
                if 'rel="next"' in part:
                    url = part.split(';')[0].strip(' <>')
                    break
    return all_results
```

---

### 3.3 RCSB PDB Search API v2 — 搜索 PDB 结构

**端点**: `POST https://search.rcsb.org/rcsbsearch/v2/query`

**Content-Type**: `application/json`

**Query 结构** — 组合 taxonomy_lineage + full_text 搜索:

```json
{
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
          "value": "10242"
        }
      },
      {
        "type": "terminal",
        "service": "full_text",
        "parameters": {
          "value": "antibody"
        }
      }
    ]
  },
  "return_type": "entry",
  "request_options": {
    "paginate": {"start": 0, "rows": 1000}
  }
}
```

**参数说明**:
- `rcsb_entity_source_organism.taxonomy_lineage.id` + `exact_match`: 匹配该 taxonomy 所有后代的 PDB (已验证匹配 10245 等后代)
- `full_text` + `value`: 全文搜索，**空格分隔 = OR 逻辑** (contains_words, any-of match)
- `rows`: 单次最多 1000 条 (已验证)
- 对于 `add_keywords` 多个关键词，拼接为空格分隔的字符串: `"antibody antigen fab"`

**返回结构**:
```json
{
  "total_count": 36,
  "result_set": [{"identifier": "2I9L", "score": 0.98}, ...]
}
```

**分页**: 若 `total_count > rows`，调整 `request_options.paginate.start` 循环直至取完。

---

### 3.4 RCSB PDB FASTA 下载

**端点**: `GET https://www.rcsb.org/fasta/entry/{PDB_ID}/download`

**返回**: 多序列 FASTA 文件，每个 chain 一个条目。

**FASTA header 格式**:
```
>{PDB_ID}_{entity_id}|Chain(s) {chain_ids}|{protein_name}|{organism} ({tax_id})
```

示例:
```
>2I9L_1|Chains A, C, E, G|Antibody 7D11 light chain|Mus musculus (10090)
>2I9L_2|Chains B, D, F, H|Antibody 7D11 heavy chain|Mus musculus (10090)
>2I9L_3|Chains I, J, K, L|Virion membrane protein M25|Vaccinia virus (10245)
```

**批量下载策略**:
- RCSB **没有** FASTA 批量 POST 端点 (返回 404)
- 需逐个 PDB ID 下载，每个约 0.3-1s
- 建议使用 `concurrent.futures.ThreadPoolExecutor` 并发 (max_workers=5) 加速，注意远端 API 与本地网络吞吐限制

**并发下载伪代码**:
```python
from concurrent.futures import ThreadPoolExecutor, as_completed

def fetch_fasta(pdb_id):
    url = f"https://www.rcsb.org/fasta/entry/{pdb_id}/download"
    resp = requests.get(url, timeout=30)
    resp.raise_for_status()
    return pdb_id, resp.text

with ThreadPoolExecutor(max_workers=5) as executor:
    futures = {executor.submit(fetch_fasta, pid): pid for pid in pdb_ids}
    for future in as_completed(futures):
        pdb_id, fasta_text = future.result()
        # parse fasta...
```

---

### 3.5 FASTA 解析 (BioPython SeqIO)

```python
from Bio import SeqIO
from io import StringIO

def parse_pdb_fasta(fasta_text):
    """返回 [(chain_desc, sequence), ...]"""
    records = []
    for record in SeqIO.parse(StringIO(fasta_text), "fasta"):
        # record.id = "2I9L_1"
        # record.description = "2I9L_1|Chains A, C, E, G|Antibody 7D11 light chain|Mus musculus (10090)"
        records.append((record.description, str(record.seq)))
    return records
```

---

## 4. 序列比对策略

### 4.1 方案选择: BioPython PairwiseAligner (局部比对)

**为什么不用 BLAST/mmseqs2**:
1. BLAST 需要额外的数据库构建步骤
2. mmseqs2 需要额外安装且命令行耦合
3. 蛋白质组规模: target ~200 条 × PDB chains ~100 条 = 约 20,000 次比对，`PairwiseAligner` 完全够用
4. `PairwiseAligner` 返回 alignment score 与 aligned sequences，便于计算 identity

**当前实现参数**:
- `mode = "local"`
- `open_gap_score = -10.0`
- `extend_gap_score = -0.5`
- 若需要使用 BLOSUM62 等替换矩阵，应显式设置 `substitution_matrix`，不要依赖默认值

**identity 计算**:
```python
from Bio.Align import PairwiseAligner

aligner = PairwiseAligner()
aligner.mode = "local"
aligner.open_gap_score = -10.0
aligner.extend_gap_score = -0.5

def align_and_identity(seq1, seq2):
    """返回 (identity_percent, score)"""
    alignments = aligner.align(seq1, seq2)
    if not alignments:
        return 0.0, 0.0
    best = alignments[0]
    aligned1 = str(best[0])
    aligned2 = str(best[1])
    score = best.score
    # 计算 identity: 相同残基数 / 比对长度
    matches = sum(1 for a, b in zip(aligned1, aligned2) if a == b)
    aligned_len = len(aligned1)  # 含 gap
    identity = matches / aligned_len if aligned_len > 0 else 0.0
    return identity * 100, score
```

**性能优化**:
- 先按序列长度过滤: 若 query 和 target 长度差 > 5x，跳过
- 只保留 identity ≥ `threshold` (默认 30%) 的结果
- 当总比对数超过阈值时，先用 k-mer 快速预筛 (Jaccard index)，仅对候选进行局部比对

### 4.2 可选: k-mer 预筛 (加速大规模比对)

```python
def kmer_jaccard(seq1, seq2, k=3):
    """3-mer Jaccard index: 快速估算序列相似性"""
    kmers1 = {seq1[i:i+k] for i in range(len(seq1)-k+1)}
    kmers2 = {seq2[i:i+k] for i in range(len(seq2)-k+1)}
    intersection = len(kmers1 & kmers2)
    union = len(kmers1 | kmers2)
    return intersection / union if union > 0 else 0.0
```

当 total_alignments > 10000 时启用预筛: Jaccard ≥ 0.1 才进入局部比对。

---

## 5. CSV 输出格式

### 5.1 Columns

| Column | 说明 | 示例 |
|--------|------|------|
| `target_taxonomy_id` | 目标物种 taxonomy ID | `10253` |
| `parent_taxonomy_id` | 父类群 taxonomy ID | `10242` |
| `target_protein_uniprot_id` | UniProt accession | `Q9JFA8` |
| `target_protein_name` | 蛋白推荐名 | `mRNA-capping enzyme catalytic subunit` |
| `support_pdb_id` | 有同源匹配的 PDB ID (分号分隔) | `2I9L;5EOR;4U6H` |
| `support_pdb_chain_description` | PDB chain 描述 (来自 FASTA header) | `2I9L_3\|Chains I,J,K,L\|Virion membrane protein M25\|Vaccinia virus (10245)` |
| `identity_pct` | 最高 identity (百分比) | `45.2` |
| `alignment_score` | 比对得分 | `823.5` |

**输出逻辑**: 每个 target protein 一行。若同一蛋白在多个 PDB 中均有匹配，`support_pdb_id` 列用分号聚合各自的最高 identity chain。

### 5.2 聚合规则

```
For each target_protein:
    matched_pdbs = []
    For each pdb_chain in all_pdb_chains:
        identity, score = align(target_seq, pdb_seq)
        if identity >= threshold:
            matched_pdbs.append((pdb_id, chain_desc, identity, score))
    # 按 identity 降序, 取每个 PDB 的最高 identity chain
    best_per_pdb = {}
    for pdb_id, chain_desc, identity, score in matched_pdbs:
        if pdb_id not in best_per_pdb or identity > best_per_pdb[pdb_id][1]:
            best_per_pdb[pdb_id] = (chain_desc, identity, score)
    # 写入一行 CSV
    row = {
        "support_pdb_id": ";".join(best_per_pdb.keys()),
        "support_pdb_chain_description": ";".join(v[0] for v in best_per_pdb.values()),
        "identity_pct": max(v[1] for v in best_per_pdb.values()),
        ...
    }
```

### 5.3 额外输出文件

1. **`{target}_{parent}_mapping.csv`** — 主结果文件
2. **`{target}_{parent}_all_alignments.csv`** — 所有通过 threshold 的比对详情 (每个 target × pdb_chain 一行)，用于深度分析
3. **`{target}_{parent}_log.txt`** — 运行日志 (API 调用数、耗时、错误)

---

## 6. 程序结构

### 6.1 文件组织

```
protein_mapper/
├── main.py              # 入口: 解析参数, 编排流程
├── config.py            # TIMEOUT, THRESHOLD 等常量
├── taxonomy.py          # NCBI lineage 验证
├── uniprot.py           # UniProt 蛋白下载
├── pdb_search.py        # RCSB PDB 搜索 + FASTA 下载
├── alignment.py         # 序列比对引擎
├── output.py            # CSV 输出
└── utils.py             # HTTP 请求封装 (端点探测 + 重试)
```

### 6.2 主流程

```python
def main(target_taxid, parent_taxid, keywords, threshold=30.0, output_dir="."):
    # 1. 验证 taxonomy lineage
    verify_lineage(target_taxid, parent_taxid)

    # 2. 下载 target 蛋白组
    target_proteins = fetch_uniprot_proteins(target_taxid)
    log(f"Downloaded {len(target_proteins)} proteins from UniProt")

    # 3. 搜索 PDB 结构
    pdb_ids = search_pdb_structures(parent_taxid, keywords)
    log(f"Found {len(pdb_ids)} PDB entries")

    # 4. 下载 PDB FASTA
    pdb_chains = download_all_pdb_fasta(pdb_ids)
    log(f"Parsed {len(pdb_chains)} PDB chains")

    # 5. 序列比对 + 输出
    results = run_alignment(target_proteins, pdb_chains, threshold)
    write_csv(results, output_dir)
    log(f"Done: {len(results)} target proteins with PDB matches")
```

### 6.3 错误处理

- **API 超时**: 重试 3 次 (指数退避: 2s, 4s, 8s)
- **个别 PDB FASTA 404**: 跳过该 PDB，记录到日志
- **空结果**: 若 PDB 搜索或比对无结果，输出空 CSV (含表头) + WARNING 日志

---

## 7. 命令行接口

```bash
python main.py \
    --target 10253 \
    --parent 10242 \
    --keywords antibody \
    --threshold 30.0 \
    --output ./results/
```

| 参数 | 必需 | 默认值 | 说明 |
|------|------|--------|------|
| `--target` | ✓ | - | Target taxonomy ID |
| `--parent` | ✓ | - | Parent taxonomy ID |
| `--keywords` | ✓ | - | 空格分隔的关键词 (AND 逻辑) |
| `--threshold` | ✗ | 30.0 | Identity 阈值 (%) |
| `--output` | ✗ | `./` | 输出目录 |

---

## 8. API 调用汇总

| Step | API | Method | 调用次数 | 备注 |
|------|-----|--------|----------|------|
| 1 | NCBI Datasets v2 | GET | 1 | |
| 2 | UniProt search | GET | ~1 per 500 proteins | 分页 |
| 3 | RCSB PDB Search v2 | POST | ~1 per 1000 PDBs | 分页 |
| 4 | RCSB PDB FASTA | GET | N (PDB entry 数) | 并发 5 worker |
| 5 | 本地 PairwiseAligner | 本地 | M×N | M=蛋白数, N=chain数 |

---

## 9. 已验证的数据

| 测试项 | 结果 |
|--------|------|
| NCBI lineage: 10242 在 10253 的 lineage 中 | ✓ True |
| UniProt: taxonomy 10253 有 239 个条目 | ✓ |
| RCSB: parent=10242 + keyword="antibody" → 36 PDBs | ✓ |
| RCSB: parent=10242 + keyword="antibody antigen" → 10 PDBs | ✓ |
| RCSB FASTA header 含 chain + protein name + organism | ✓ |

---

## 10. 扩展方向 (可选)

- 支持 `--alignment-mode blast` 调用本地 blastp
- 支持 `--kmer-pre filter` 切换预筛策略
- 输出 HTML 可视化报告
- 支持多个 target taxonomy 批量比对
