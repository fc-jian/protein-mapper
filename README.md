# Protein Sequence Homology Mapping Pipeline

自动化跨物种蛋白质序列同源映射工具。给定一个目标物种（target）和一个父类群（parent），从 UniProt 抓取目标全蛋白质组，从 RCSB PDB 搜索父类群下含指定关键词的结构，下载 FASTA 序列后进行全量序列比对，输出同源蛋白 ↔ PDB 结构对应关系 CSV。

## 工作流程

```
输入: target_taxid, parent_taxid, add_keywords[]
  │
  ├─ [1] NCBI Datasets v2 → 验证 parent 是否包含 target（lineage 血缘检查）
  ├─ [2] UniProt REST API → 下载 target 全蛋白质组完整序列（~239 proteins）
  ├─ [3] RCSB PDB Search v2 → 搜索 parent 下含有所有关键词的 PDB 结构
  ├─ [4] RCSB PDB FASTA → 并发下载每个 PDB 的链序列
  └─ [5] PairwiseAligner + k-mer 预筛 → 全量比对，输出 CSV
```

## 快速开始

### 环境要求

- Python 3.10+
- 可访问 NCBI、UniProt 与 RCSB PDB 的网络环境

### 安装

```bash
cd /path/to/protein_mapper
python3 -m venv .venv
source .venv/bin/activate
pip install pandas requests biopython
```

### 运行

```bash
source .venv/bin/activate
python main.py \
    --target 10253 \           # 目标物种 taxonomy ID（Vaccinia virus Tian Tan）
    --parent 10242 \           # 父类群 taxonomy ID（Chordopoxvirinae）
    --keywords antibody \      # PDB 搜索关键词（AND 逻辑，可多个）
    --threshold 30.0 \         # Identity 阈值 (%)
    --output ./results/        # 输出目录
```

## 命令行参数

| 参数 | 必需 | 默认值 | 说明 |
|------|:----:|--------|------|
| `--target` | ✓ | — | Target taxonomy ID |
| `--parent` | ✓ | — | Parent taxonomy ID |
| `--keywords` | ✓ | — | 空格分隔的关键词，AND 逻辑 |
| `--threshold` | | 30.0 | Identity 百分比阈值 |
| `--output` | | `./` | 输出目录 |

## 输出文件

| 文件 | 说明 |
|------|------|
| `{target}_{parent}_mapping.csv` | 核心映射表（每个 target protein 一行） |
| `{target}_{parent}_all_alignments.csv` | 全部通过阈值的比对详情 |
| `{target}_{parent}_log.txt` | 运行日志与统计摘要 |

### mapping.csv 列说明

| 列 | 示例 |
|----|------|
| `target_taxonomy_id` | `10253` |
| `parent_taxonomy_id` | `10242` |
| `target_protein_uniprot_id` | `Q9JFA1` |
| `target_protein_name` | `Cell surface-binding protein OPG105` |
| `target_organism` | `Vaccinia virus (strain Tian Tan)` |
| `support_pdb_id` | `4E9O;4ETQ;5USH;5USL;6B9J` |
| `support_pdb_chain_description` | 各 PDB chain 的 FASTA header 描述 |
| `identity_pct` | `83.28` |
| `alignment_score` | `259.0` |

当同一蛋白在多个 PDB 中均有匹配时，`support_pdb_id` 和 `support_pdb_chain_description` 用分号聚合各 PDB 最高 identity 的 chain。

## 示例结果

```bash
$ python main.py --target 10253 --parent 10242 --keywords antibody

[Step 2/5] 验证 NCBI Taxonomy 血缘关系
  ✓ Lineage 验证通过: 10242 是 10253 的祖先

[Step 3/5] 下载目标物种蛋白质组 (UniProt)
  ✓ UniProt 下载完成: 共 239 个蛋白质

[Step 4/5] 搜索 RCSB PDB 结构
  ✓ PDB 搜索完成: 共 36 个结构

[Step 5/5] 序列比对
  开始序列比对: 239 目标蛋白 × 114 PDB chains = 27246 次比对
  启用 k-mer (k=3) 预筛, Jaccard ≥ 0.1
  比对统计: 长度过滤 1766, k-mer 过滤 25431, 匹配 49
  ✓ 比对完成: 10 个目标蛋白有 PDB 匹配

✓ 管道执行完成! 耗时 18.2s
  统计: 239 目标蛋白, 36 PDB, 10 蛋白有匹配
```

### Top 匹配

| Identity | Target Protein | PDB Support |
|----------|---------------|-------------|
| **83.3%** | Cell surface-binding protein OPG105 | 4E9O, 4ETQ, 5USH, 5USL, 6B9J |
| **81.4%** | TK3L | 9HBK, 9HL2, 9RDH |
| **76.5%** | Protein OPG190 | 8XS3 |
| **76.2%** | Virion membrane protein A16 | 9HBK~9RDH (8 PDBs) |
| **75.9%** | Entry-fusion complex protein OPG094 | 9HBK~9RDH (8 PDBs) |

## API 端点一览

| 步骤 | API | Method | 端点 |
|------|-----|--------|------|
| 1 | NCBI Datasets v2 | GET | `api.ncbi.nlm.nih.gov/datasets/v2/taxonomy/taxon/{taxid}` |
| 2 | UniProt REST | GET | `rest.uniprot.org/uniprotkb/search?query=taxonomy_id:{taxid}` |
| 3 | RCSB PDB Search v2 | POST | `search.rcsb.org/rcsbsearch/v2/query` — `contains_words` (any-of / OR) |
| 4 | RCSB PDB FASTA | GET | `www.rcsb.org/fasta/entry/{pdb_id}/download` |

## 文件结构

```
protein_mapper/
├── main.py            # 入口：解析参数，编排流程
├── config.py          # 超时、并发、阈值等常量
├── utils.py           # HTTP 请求封装（端点探测 + 指数退避重试）
├── taxonomy.py        # NCBI lineage 验证
├── uniprot.py         # UniProt 蛋白质组下载（分页处理）
├── pdb_search.py      # RCSB PDB 搜索 + 并发 FASTA 下载
├── alignment.py       # PairwiseAligner + k-mer 预筛引擎
├── output.py          # CSV 输出
├── docs/
│   └── DEVELOPMENT_DOC.md # 完整技术文档
└── README.md          # 本文件
```

## 配置常量

所有可调参数在 `config.py` 中：

```python
TIMEOUT = 60                # 常规请求超时 (s)
FASTA_TIMEOUT = 30          # PDB FASTA 下载超时 (s)
MAX_RETRIES = 3             # 重试次数（指数退避: 2s→4s→8s）
UNIPROT_PAGE_SIZE = 500     # UniProt 分页大小
PDB_PAGE_SIZE = 1000        # PDB 搜索分页大小
PDB_CONCURRENT_WORKERS = 5  # PDB FASTA 并发下载数
DEFAULT_THRESHOLD = 30.0    # 默认 identity 阈值 (%)
LENGTH_RATIO_MAX = 5.0      # 序列长度比过滤上限
KMER_SIZE = 3               # k-mer 预筛 k 值
KMER_JACCARD_THRESHOLD = 0.1  # Jaccard 预筛阈值
KMER_PRE_FILTER_TRIGGER = 10000  # 触发 k-mer 预筛的比对总数
```

## 依赖

| 库 | 用途 |
|----|------|
| `requests` | HTTP 请求 |
| `biopython` | PairwiseAligner + SeqIO（FASTA 解析） |
| `pandas` | CSV 读写 |

## 许可

内部研究工具。完整技术文档见 `DEVELOPMENT_DOC.md`。
