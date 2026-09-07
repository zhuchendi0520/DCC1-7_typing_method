# DCC1-DCC7 typing for *Mycobacterium abscessus*
# 脓肿分枝杆菌 DCC1-DCC7 分型

This repository provides an evolutionary path-based SNP barcode workflow for assigning *Mycobacterium abscessus* isolates to the seven recognized dominant circulating clones (DCC1-DCC7). It contains the two required reference genomes, curated DCC-specific marker sets, scripts for converting VarScan consensus output, a unified assignment program, and an optional end-to-end FASTQ workflow.

本方法包基于进化路径上的 SNP 标记，对脓肿分枝杆菌进行 DCC1-DCC7 分型，包含两套参考基因组、筛选后的 DCC 特异性标记、VarScan 结果转换脚本、统一分型程序及从双端 FASTQ 开始的一体化流程。

## Scope / 适用范围

The typing scheme supports:

| Subspecies / 亚种 | Reference genome / 参考基因组 | DCCs evaluated / 分型范围 |
|---|---|---|
| subsp. *abscessus* (ABS) | GZ002, NZ_CP034181.1 | DCC1, DCC2, DCC4, DCC5 |
| subsp. *massiliense* (MAS) | CCUG 48898, NZ_AP014547.1 | DCC3, DCC6, DCC7 |

Subsp. *bolletii* is not included because DCC1-DCC7 do not belong to this subspecies. 

本工具仅针对 DCC1-DCC7，不包含 *bolletii* 亚种。

## How the marker sets were derived / 标记集的来源

The following summarizes the derivation described in the accompanying study; the packaged scripts apply the resulting marker sets rather than repeat marker discovery.

以下概述配套研究中的标记构建过程。本包脚本用于应用已经得到的标记集，不重新执行标记发现流程。

1. Subspecies-specific phylogenies were reconstructed from whole-genome sequence data after masking recombinant regions.
   利用全基因组序列，在屏蔽重组区域后分别重建各亚种的系统发育树。
2. The seven previously recognized DCC clades were located on the corresponding maximum-likelihood trees.
   在对应的最大似然树上定位已知的七个 DCC 分支。
3. Large DCCs were pruned with Treemmer while retaining at least 70% of their phylogenetic diversity.
   对较大的 DCC 使用 Treemmer 缩减样本，同时保留至少 70% 的系统发育多样性。
4. SNPPar v1.0 was used to reconstruct ancestral states under a GTR model with gamma-distributed rate heterogeneity.
   使用 SNPPar v1.0 重建祖先状态，采用 GTR 模型及位点间速率异质性的 gamma 分布。
5. Mutations accumulated along the evolutionary path leading to each DCC ancestor were collected as candidate barcode SNPs.
   提取通向各 DCC 祖先的进化路径上累积的突变，作为候选 SNP 条形码。
6. Candidate SNPs shared nonspecifically across lineages were removed.
   去除不同谱系间非特异性共享的候选 SNP。
7. A marker was retained when its derived allele was present in at least 95% of isolates within the target DCC and in less than 5% of isolates outside that DCC.
   保留目标 DCC 内携带率至少为 95%、目标 DCC 外携带率低于 5% 的衍生等位基因标记。

The final marker counts are / 最终标记数量如下：

| DCC | Subspecies / 亚种 | Markers / 标记数 |
|---|---|---:|
| DCC1 | ABS | 819 |
| DCC2 | ABS | 578 |
| DCC3 | MAS | 257 |
| DCC4 | ABS | 478 |
| DCC5 | ABS | 154 |
| DCC6 | MAS | 602 |
| DCC7 | MAS | 1,184 |

The full article-ready description is provided in [`docs/METHODS.md`](docs/METHODS.md).

完整英文方法说明见 [`docs/METHODS.md`](docs/METHODS.md)。

## Directory structure / 文件结构

```text
DCC1-7_typing_method/
├── README.md
├── docs/
│   └── METHODS.md
├── examples/
│   └── expected_output.tsv
├── markers/
│   ├── ABS_DCC_markers.tsv
│   └── MAS_DCC_markers.tsv
├── references/
│   ├── ABS_GZ002_NZ_CP034181.1.fna
│   ├── MAS_CCUG48898_NZ_AP014547.1.fna
│   └── SHA256SUMS
└── scripts/
    ├── 01_cns_to_typesnp.py
    ├── 02_assign_dcc.py
    └── 03_run_dcc_typing.sh
```

`markers/` 保存分型标记，`references/` 保存对应参考序列及校验值，`scripts/` 保存编号脚本，`examples/` 提供输出格式示例。

## Requirements / 运行环境

- Python 3.9 or later
- BWA 0.7.17 or later
- SAMtools
- Java
- VarScan 2.3.9
- Paired-end Illumina FASTQ files for the end-to-end workflow

The two Python scripts use only the Python standard library.

两个 Python 脚本仅使用标准库，需要 Python 3.9 或以上版本。从 FASTQ 开始的完整流程还需要 BWA、SAMtools、Java 和 VarScan，以及双端 Illumina 测序数据。请先完成测序数据质量控制。

## Important coordinate-system requirement / 参考坐标必须一致

The sample calls and marker set must use the same reference coordinate system. Map ABS isolates to GZ002 and MAS isolates to CCUG 48898. Calls generated against a different reference cannot be compared directly with these marker files.

The isolate's subspecies must therefore be determined before DCC assignment, for example by average nucleotide identity against representative genomes.

**样本与标记集必须使用相同的参考坐标。** ABS 样本使用 GZ002，MAS 样本使用 CCUG 48898。不能直接将基于 ATCC19977 或其他参考序列得到的位点用于本包标记集。请在分型前完成亚种鉴定，例如通过与代表性基因组计算 ANI 确定亚种。

## Option 1: End-to-end typing from paired FASTQ files / 方案一：从双端 FASTQ 开始

Run the commands below from the package root and replace the example paths with your own paths.

以下命令均在本方法包根目录执行，请将示例路径替换为实际文件路径。

Make the scripts executable / 添加执行权限：

```bash
chmod +x scripts/*.py scripts/*.sh
```

Run an ABS isolate / 对 ABS 样本运行：

```bash
scripts/03_run_dcc_typing.sh \
  --sample sample01 \
  --r1 /path/to/sample01_R1.fastq.gz \
  --r2 /path/to/sample01_R2.fastq.gz \
  --subspecies ABS \
  --varscan-jar /path/to/VarScan.v2.3.9.jar \
  --threads 8 \
  --outdir results/sample01
```

For a MAS isolate, use `--subspecies MAS`. The final result is written to:

MAS 样本将参数改为 `--subspecies MAS`，程序会选择对应参考序列和标记集。最终结果输出至：

```text
results/sample01/sample01.dcc_assignment.tsv
```

The workflow performs the following operations:

1. BWA-MEM mapping to the appropriate reference.
2. BAM sorting and indexing with SAMtools.
3. Pileup generation with mapping quality >=30 and base quality >=20.
4. VarScan consensus calling using a minimum variant frequency of 0.75.
5. Conversion to a position/allele table.
6. DCC scoring and assignment.

流程依次执行：BWA-MEM 比对、SAMtools 排序和索引、生成 pileup（比对质量至少 30，碱基质量至少 20）、VarScan 共识调用（最低变异频率 0.75）、转换位点/等位基因表，以及 DCC 评分与分型。

## Option 2: Type an existing VarScan CNS file / 方案二：使用已有 CNS 文件

Convert the CNS file / 转换 VarScan CNS 文件：

```bash
python3 scripts/01_cns_to_typesnp.py sample01.cns \
  --output sample01.typesnp.tsv
```

Assign an ABS isolate / 对 ABS 样本分型：

```bash
python3 scripts/02_assign_dcc.py \
  markers/ABS_DCC_markers.tsv \
  sample01.typesnp.tsv \
  --subspecies ABS \
  --sample sample01 \
  > sample01.dcc_assignment.tsv
```

Assign a MAS isolate by replacing the marker file and subspecies:

MAS 样本需要同时替换标记文件和亚种参数：

```bash
python3 scripts/02_assign_dcc.py \
  markers/MAS_DCC_markers.tsv \
  sample01.typesnp.tsv \
  --subspecies MAS \
  --sample sample01 \
  > sample01.dcc_assignment.tsv
```

## Assignment rule / 分型规则

For each eligible DCC, the program calculates:

```text
detected lineage-defining derived alleles / total lineage markers x 100
```

The lineage with the highest score is reported as the best match. A sample is assigned to that DCC when the score is at least 95%; otherwise, its assignment is `Non_DCC`. Exact ties are retained as separate output rows and should be reviewed manually.

对对应亚种的每个 DCC，计算“匹配的衍生等位基因标记数 / 该 DCC 标记总数 × 100”。最高分达到 95% 时，归入对应 DCC；否则输出 `Non_DCC`。若最高分完全相同，则保留多行并列结果，需要人工复核。

Missing calls are counted as unmatched markers; the denominator is the full marker set, not only callable sites. The score is a marker-match percentage, not an assignment probability.

缺失或无法判定的位点按未匹配处理，分母始终为完整标记集的数量，而非成功调用的位点数。匹配百分比不等同于分型概率。

## Output columns / 输出字段

| Column / 字段 | Description / 含义 |
|---|---|
| `Sample` | User-supplied or inferred sample identifier / 用户提供或程序推定的样本编号 |
| `Subspecies` | ABS or MAS / 亚种 |
| `Best_match` | DCC with the highest marker match percentage / 匹配率最高的 DCC |
| `Match_percentage` | Percentage of markers carrying the expected derived allele / 与预期衍生等位基因匹配的标记百分比 |
| `Matched_markers` | Number of matched markers / 匹配标记数量 |
| `Total_markers` | Number of markers in the best-matching DCC barcode / 最佳匹配 DCC 的标记总数 |
| `Assignment` | DCC1-DCC7 if the score is >=95%; otherwise `Non_DCC` / 达到阈值则为对应 DCC，否则为 `Non_DCC` |

Example / 输出示例：

```text
Sample    Subspecies    Best_match    Match_percentage    Matched_markers    Total_markers    Assignment
sample01  ABS           DCC1          100.00              819                819              DCC1
```

## Batch processing / 批量处理

Run `03_run_dcc_typing.sh` once per isolate after separating samples by subspecies. Assignment tables can be combined after retaining one header:

按亚种区分样本后，对每个菌株分别运行 `03_run_dcc_typing.sh`。以下命令合并结果，并仅保留一行表头：

```bash
awk 'FNR == 1 && NR != 1 {next} {print}' results/*/*.dcc_assignment.tsv \
  > DCC1-7_assignments.tsv
```

## Interpretation and quality control / 结果解释与质量控制

- `Non_DCC` means that none of DCC1-DCC7 reached the 95% threshold; it does not identify a novel lineage.
- A low score can result from poor coverage, contamination, an incorrect subspecies/reference choice, or genuine non-DCC ancestry.
- Review mean depth, genome coverage, and contamination before interpreting borderline assignments. The current output reports only the best match and exact ties, not all DCC scores.
- This method is intended for lineage assignment and surveillance, not for inferring recent transmission between isolates.
- The barcode was derived from the population analyzed in the accompanying study and should be revalidated as geographically and temporally broader datasets become available.

- `Non_DCC` 仅表示没有达到这七个 DCC 的分型阈值，不代表发现新谱系。
- 低匹配率可能源于覆盖不足、污染、亚种或参考序列选错，也可能确实为非 DCC 菌株。
- 对接近阈值的样本，应检查平均深度、基因组覆盖率和污染情况。当前结果仅报告最佳匹配及完全并列的匹配，不输出所有 DCC 的评分。
- 本工具用于谱系分型和监测，不能仅据此推断近期传播关系。
- 标记来源于配套研究所覆盖的群体，随着地域和时间范围扩大，应在新增数据中继续验证。

## Citation / 引用

When using this workflow, cite the accompanying study and the software used for mapping, consensus calling, phylogenetic reconstruction, ancestral reconstruction, and diversity-preserving subsampling.

使用本流程时，请引用配套研究，以及实际使用的比对、共识调用、系统发育重建、祖先重建和保留多样性抽样软件。正式发布时，请补充配套研究的完整引用信息。
