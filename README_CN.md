# 蛋白质序列同源映射管道

[English](README.md)

Protein Mapper 是一个 Python 命令行工具，用于自动化完成跨物种蛋白质序列同源映射。给定目标物种、父级分类群 ID 和 PDB 搜索关键词后，程序会从 UniProt 获取目标全蛋白质组，从 RCSB PDB 搜索父级分类群下的相关结构，下载 PDB FASTA 链序列，执行序列比对，并输出 CSV 映射结果。

## 工作流程

```text
输入: target_taxid, parent_taxid, keywords[]
  |
  |-- [1] NCBI Datasets v2: 验证 parent lineage 是否包含 target
  |-- [2] UniProt REST API: 获取目标物种蛋白质组
  |-- [3] RCSB PDB Search API v2: 按父级 lineage 和关键词搜索结构
  |-- [4] RCSB PDB FASTA: 下载并解析 PDB chain 序列
  |-- [5] PairwiseAligner + k-mer 预筛: 序列比对
  `-- [6] MMseqs2: 对 mapped UniProt sequences 聚类并输出报告
```

## 环境要求

- Python 3.10+
- `uv`
- MMseqs2（`mmseqs` 在 `PATH` 中可用）
- 可访问 NCBI、UniProt 和 RCSB PDB 公开 API 的网络环境

## 安装

本仓库使用 `uv` 管理本地环境，虚拟环境位于 `./.venv`。

```bash
cd /path/to/protein-mapper
uv venv .venv --python 3.10
source .venv/bin/activate
uv pip install requests biopython pandas
```

MMseqs2 需要通过系统或生信软件包管理器单独安装，并确认 `mmseqs` 命令可用。

## 运行

```bash
source .venv/bin/activate
python main.py \
    --target 10253 \
    --parent 10242 \
    --keywords antibody Fab scFv nanobody VHH \
    --threshold 30.0 \
    --cluster-threshold 0.9 \
    --output ./results/
```

## 命令行参数

| 参数 | 必需 | 默认值 | 说明 |
| --- | :---: | --- | --- |
| `--target` | 是 | - | 目标 NCBI taxonomy ID |
| `--parent` | 是 | - | 用于 lineage 过滤的父级 NCBI taxonomy ID |
| `--keywords` | 是 | - | 空格分隔的 PDB 搜索关键词 |
| `--threshold` | 否 | `30.0` | 最低 identity 百分比 |
| `--cluster-threshold` | 否 | `0.9` | 对 mapped UniProt proteins 运行 MMseqs2 聚类时使用的序列 identity 阈值 |
| `--output` | 否 | `./results/` | 输出目录 |

## 输出文件

每次运行会生成四个报告文件：

| 文件 | 说明 |
| --- | --- |
| `{target}_{parent}_mapping.csv` | 主映射表，每个匹配到 PDB 的 target protein 一行 |
| `{target}_{parent}_all_alignments.csv` | 所有通过 identity 阈值的比对详情 |
| `{target}_{parent}_cluster_mapping.csv` | 对 mapped UniProt sequences 进行 MMseqs2 聚类后的 cluster 级映射表 |
| `{target}_{parent}_log.txt` | 本次运行的统计摘要 |

`results/` 属于生成结果目录，不由 Git 跟踪。

程序默认保存中间 FASTA 序列和元数据，位置如下：

```text
{output}/{target}_{parent}_intermediates/
├── uniprot/
│   ├── {target}_uniprot.fasta
│   └── {target}_uniprot_proteins.csv
├── pdb/
│   ├── {target}_{parent}_pdb_ids.txt
│   ├── {target}_{parent}_pdb_query.json
│   ├── {target}_{parent}_pdb_chains.csv
│   └── fasta/{pdb_id}.fasta
└── mmseqs/
    └── {target}_{parent}_mapped_uniprot.fasta
```

Cluster report 会按配置的序列 identity 阈值聚合最终 mapping 通过的 UniProt sequences。每个 cluster 内，所有 mapping columns 都会先去重，再用分号连接。

### Mapping 列说明

| 列 | 示例 |
| --- | --- |
| `target_taxonomy_id` | `10253` |
| `parent_taxonomy_id` | `10242` |
| `target_protein_uniprot_id` | `Q9JFA1` |
| `target_protein_name` | `Cell surface-binding protein OPG105` |
| `target_organism` | `Vaccinia virus (strain Tian Tan)` |
| `support_pdb_id` | `4E9O;4ETQ;5USH;5USL;6B9J` |
| `support_pdb_chain_description` | PDB FASTA chain 描述 |
| `identity_pct` | `83.28` |
| `alignment_score` | `259.0` |

当同一个 target protein 匹配多个 PDB 条目时，`support_pdb_id` 和 `support_pdb_chain_description` 使用分号连接。聚合结果中，每个 PDB ID 保留 identity 最高的 chain。

## 示例

```bash
python main.py --target 10253 --parent 10242 --keywords antibody --output ./results/
```

成功运行时会按步骤输出摘要：

```text
[Step 1/6] Verify NCBI Taxonomy lineage
[Step 2/6] Download target proteome from UniProt
[Step 3/6] Search RCSB PDB structures
[Step 4/6] Download PDB FASTA sequences
[Step 5/6] Sequence alignment
[Step 6/6] Cluster mapped UniProt sequences with MMseqs2

Pipeline completed.
```

## API 端点

| 步骤 | API | Method | Endpoint |
| --- | --- | --- | --- |
| 1 | NCBI Datasets v2 | GET | `api.ncbi.nlm.nih.gov/datasets/v2/taxonomy/taxon/{taxid}` |
| 2 | UniProt REST | GET | `rest.uniprot.org/uniprotkb/search?query=taxonomy_id:{taxid}` |
| 3 | RCSB PDB Search v2 | POST | `search.rcsb.org/rcsbsearch/v2/query` |
| 4 | RCSB PDB FASTA | GET | `www.rcsb.org/fasta/entry/{pdb_id}/download` |

## 仓库结构

```text
protein-mapper/
├── main.py        # CLI 入口和流程编排
├── config.py      # 超时、分页、并发和阈值常量
├── utils.py       # HTTP 端点探测、重试和日志
├── taxonomy.py    # NCBI lineage 验证
├── uniprot.py     # UniProt 蛋白质组下载和 JSON 解析
├── pdb_search.py  # RCSB PDB 搜索和 FASTA 下载/解析
├── alignment.py   # PairwiseAligner + k-mer 预筛
├── intermediate.py # 中间 FASTA 和元数据输出
├── clustering.py  # MMseqs2 聚类和 cluster 级报告生成
├── output.py      # CSV 和摘要日志输出
├── AGENTS.md      # 协作者和代理说明
├── README.md      # 英文文档
└── README_CN.md   # 中文文档
```

## 配置

常用可调参数位于 `config.py`：

```python
TIMEOUT = 60
FASTA_TIMEOUT = 30
MAX_RETRIES = 3
SAVE_INTERMEDIATES = True
INTERMEDIATE_DIR_SUFFIX = "_intermediates"
UNIPROT_PAGE_SIZE = 500
PDB_PAGE_SIZE = 1000
PDB_CONCURRENT_WORKERS = 5
DEFAULT_THRESHOLD = 30.0
LENGTH_RATIO_MAX = 5.0
KMER_SIZE = 3
KMER_JACCARD_THRESHOLD = 0.1
KMER_PRE_FILTER_TRIGGER = 10000
MMSEQS_BINARY = "mmseqs"
MMSEQS_CLUSTER_THRESHOLD = 0.9
MMSEQS_THREADS = 1
```

## 依赖

| 包 | 用途 |
| --- | --- |
| `requests` | HTTP 请求 |
| `biopython` | 序列比对和 FASTA 解析 |
| `pandas` | CSV 输出 |
| MMseqs2 | 对 mapped UniProt sequences 聚类 |

## 许可

内部研究工具。
