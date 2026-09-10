# `--stream-by-region` 与 GATK / 非流式输出不一致（2026-09-10）

状态：**已复现，P1（流式路径正确性问题）**。发现于第 1 轮 Track A 审计。
本报告只记录证据与假设，不包含修复。

## 结论摘要

`fastgatk-hc-call --stream-by-region N`（当 tile 小于区间长度时）的输出与
**同一条命令去掉 `--stream-by-region`** 的输出不一致，也与 pinned GATK 4.6.2.0 不一致。
不一致有两类：

1. **参考块粒度改变（结构性）**：非流式路径与 GATK 都把 20:10020230–10020428 拆成
   29 个小 `<NON_REF>` 块；流式路径把它合并成 **1 个** 块（`POS=10020230 END=10020428`）。
   在 tile=100/500 时，流式输出比非流式输出**少 28 条记录**。
2. **注解读数错误（数值性）**：在 20:10020680 这一行，流式路径给出的
   `RAW_MQandDP` 与 `SB` 与非流式路径、与 GATK 都不同。

第 2 类可由单行证据直接判定：非流式路径与 GATK **逐字段相同**，流式路径不同。

## 最小复现

参考/输入夹具：`fixtures/chr20/ref20mnp.fasta`、`fixtures/chr20/mnp.bam`（仓库内已有）。

```bash
B=fastgatk-native/build/fastgatk-hc-call
COMMON=(-R fixtures/chr20/ref20mnp.fasta -I fixtures/chr20/mnp.bam \
        -L 20:10019901-10020710 --emit-ref-confidence GVCF --max-mnp-distance 1 \
        --threads 1 --add-output-vcf-command-line false)

$B "${COMMON[@]}"                          -O whole.g.vcf    # 48 条记录
$B "${COMMON[@]}" --stream-by-region 300   -O t300.g.vcf     # 48 条，但字段不同
$B "${COMMON[@]}" --stream-by-region 500   -O t500.g.vcf     # 20 条
```

## 证据 1：tile 尺寸扫描（OpenMP 与 Serial 结果一致）

| `--stream-by-region` | 记录数 | 相对非流式缺失 | 相对非流式多出 | 与非流式逐字节相同 |
| --- | ---: | ---: | ---: | --- |
| 100 | 20 | 28 | 0 | 否 |
| 200 | 41 | 8 | 1 | 否 |
| 300 | 48 | 0 | 0 | **否** |
| 405 | 41 | 8 | 1 | 否 |
| 500 | 20 | 28 | 0 | 否 |
| 810（=区间长度，单 tile） | 48 | 0 | 0 | **是** |
| 1000 | 48 | 0 | 0 | **是** |
| 1000000 | 48 | 0 | 0 | **是** |

要点：**只有 tile 覆盖整个区间时，流式输出才与非流式输出逐字节相同**。tile 变小时偏离出现，
且偏离程度随 tile 尺寸非单调变化。`--threads 1` 固定，排除线程数干扰。

## 证据 2：20:10020680 处注解读数

非流式（OpenMP 与 Serial 相同）：

```text
20  10020680  .  CA  AT,*,<NON_REF>  154.25  .  DP=8;ExcessHet=0.0000;MLEAC=1,1,0;MLEAF=0.500,0.500,0.00;RAW_MQandDP=28800,8  GT:AD:DP:GQ:PGT:PID:PL:PS:SB  1|2:0,4,2,0:6:84:1|0:10020679_AC_TA:249,84,99,165,0,159,252,99,168,264:10020679:0,0,3,3
```

流式 `--stream-by-region 300`：

```text
20  10020680  .  CA  AT,*,<NON_REF>  154.25  .  DP=8;ExcessHet=0.0000;MLEAC=1,1,0;MLEAF=0.500,0.500,0.00;RAW_MQandDP=97200,27  GT:AD:DP:GQ:PGT:PID:PL:PS:SB  1|2:0,4,2,0:6:84:1|0:10020679_AC_TA:249,84,99,165,0,159,252,99,168,264:10020679:0,0,0,0
```

差异字段：`RAW_MQandDP`（`28800,8` → `97200,27`）与 `SB`（`0,0,3,3` → `0,0,0,0`）。
`QUAL`、`PL`、`GT`、`GQ`、`PGT/PID/PS` 完全相同，即**基因型数值未受影响，受影响的是注释**。

## 证据 3：GATK 4.6.2.0 基准

```bash
third_party/jdk17/bin/java -Xmx1g -jar \
  third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar HaplotypeCaller \
  -R fixtures/chr20/ref20mnp.fasta -I fixtures/chr20/mnp.bam -L 20:10019901-10020710 \
  -ERC GVCF --max-mnp-distance 1 -O gatk.g.vcf --native-pair-hmm-threads 1 \
  --create-output-variant-index false --add-output-vcf-command-line false
```

- GATK 记录数 **48**，与非流式 native 相同。
- GATK 在 20:10020680 的行为
  `... RAW_MQandDP=28800,8 ... SB=0,0,3,3`，**与非流式 native 逐字段相同**。

结论方向明确：**非流式路径是对的（与 GATK 一致），`--stream-by-region` 是偏离的一方。**

## 与现有绿色测试的关系

`fastgatk-native/scripts/verify_hc_region_streaming.py`（CTest
`fastgatk-hc-region-streaming-contract`）断言：

- `--stream-by-region 1000000` 与普通输出**逐字节相同**；
- `--stream-by-region 500` 与普通输出**locus key 集合相同**。

该测试目前在双后端全绿，但它使用的是另一份夹具/区间；在
`fixtures/chr20/mnp.bam` + `20:10019901-10020710` 上，tile=500 的 key 集合
比非流式少 28 条。也就是说：**该契约在当前实现上并非普遍成立，而是依赖夹具**。
本次发现的用例应补进该 oracle，作为回归守卫。

## 根因分析（第 2 轮进展）

### 已排除：`RegionGvcfStitcher::can_merge` 过宽

初版假设是跨 tile 拼接器合并了本应独立的参考块。**该假设已被数值证据证伪**：
`RegionGvcfStitcher::merge()`（`fastgatk-native/src/hc_call.cpp:5864`）对 PL 取逐元素最小值。
若 tile=500 的那一条 28 块合并记录真是拼接器产生的，其值应为
`PL=0,26,494`、`GQ=26`、`MIN_DP=12`（对 29 个块的逐元素最小值）。

实测 tile=500 该记录为：

```text
10020230  ...  END=10020428  ...  GT:DP:GQ:MIN_DP:PL   0/0:26:0:21:0,0,0
```

即 `PL=0,0,0`、`GQ=0`、`MIN_DP=21` —— 与拼接器合并的签名不符。**拼接器不是本例的成因**，
虽其 `can_merge` 只比较 `pl.size()` 而不比较 PL 值这一点仍值得单独复核（见「未证明」）。

### 已确认：退化的参考置信度来自 tile 自身

`PL=0,0,0 / GQ=0` 意味着这些 locus **没有参考置信度似然值**，而不是"被合并"。
配合：`--stream-by-region 300` 时该窗口完全落在单个 tile 内部，输出与 GATK 的 29 块结构
**逐条相同**；只有 tile 边界切过该窗口时才退化。因此成因是
**tile 边界处的参考置信度构造**，不是跨 tile 合并。

### 候选机制（尚未证实到具体代码路径）

参考块有两个构造器：

| 构造器 | 位置 | 输入 | 被谁调用 |
| --- | --- | --- | --- |
| `build_reference_blocks`（精确） | `calling_pipeline.cpp:3295` | `rcm_loci` + `reference_confidence_observations` + `corrected_reads`（locus 级） | `calling_pipeline.cpp:15214`（`assembly_region_independent_pass` 分支） |
| `build_profile_local_reference_blocks`（近似） | `calling_pipeline.cpp:13111` | `activity.profile_regions` + `active_reads`（分区级） | `calling_pipeline.cpp:15206`（无 calling region 时）与 `15672`（partition-merge 路径） |

近似构造器使用**分区级**的 profile regions 与 `active_reads`，若它按分区而非按 locus 推导
参考置信度，则正好同时解释两类症状：

- 一个 profile region 覆盖整段 → 单块、`PL=0,0,0`；
- 注释反映分区读取集合 → `RAW_MQandDP=97200,27`（分区 27 条）与 `SB=0,0,0,0`。

**尚未证明**流式路径实际走到上述哪一个调用点。已尝试用调试开关定位：
`FASTGATK_DEBUG_GVCF_EMISSION=1` 会打印 `[FASTGATK_GVCF_BLOCK_INPUT] stage=...`
与 `[FASTGATK_GVCF_EMISSION] pos=... symbolic-pl=unavailable`；在本夹具上
**流式与非流式都只出现 `stage=direct`**，未触发 `stage=merged`，因此该实验对定位无结论。
注意 `symbolic-pl=unavailable` 与观察到的 `PL=0,0,0` 可能同源，值得优先追。

### 与状态文档旧记录的冲突

`IMPLEMENTATION_STATUS.md` 曾记录「流式 gVCF 输出路径已将实际 read batch 传给注释重算，
而非使用空 read 指针」，并据此视为已修正。本次结果表明：该改动**可能传入了非空但范围错误的
read batch（分区级而非 locus 级）**，因此"不再用空指针"并不等于"数值正确"。

## 已证明 / 未证明

**已证明**

- 非流式 native 与该区间上的 GATK 4.6.2.0 **全部 48 条数据行逐字段相同**（`diff` 无输出），
  包括 10020231–10020428 的 29 块碎片结构。这是本报告的基准：非流式路径是正确的一方。
- 在 `fixtures/chr20/mnp.bam` + `20:10019901-10020710` 上，`--stream-by-region N`
  （N < 区间长度）与非流式输出不一致：记录数最多少 28 条；tile=300 时字段
  `RAW_MQandDP`/`SB` 不同；tile=500 时出现 `PL=0,0,0 / GQ=0` 的退化块。
- 偏离只在 tile 边界切过该窗口时出现（tile=300 完全落在窗口外 → 与 GATK 逐条相同）。
- OpenMP 与 Serial 的流式偏离完全相同（Host 路径问题，非 kernel）。
- 覆盖范围本身未丢失：tile=500 时 28 条被合并进 `10020230–10020428` 单块，
  10 条变异记录逐一保留。
- **已证伪**：拼接器 `can_merge` 过宽不是本例成因（数值签名不符，见上文）。

**未证明**

- 未定位到具体代码行；`build_profile_local_reference_blocks` 与
  `build_reference_blocks` 的取舍只是候选机制，未证实流式路径走到哪一个调用点。
- `RegionGvcfStitcher::can_merge` 只比较 `pl.size()` 而不比较 PL 值：本例未触发，
  但在其他夹具上仍可能造成不正确的块合并，需单独验证。
- 未在该区间之外确定偏离边界（哪些夹具会触发、是否所有多重叠 indel 窗口都触发）。
- 未验证 `--stream-by-contig`、以及与 `--max-alternate-alleles`、read-overlap margin
  等选项的组合。
- 未验证下游影响（如 GenotypeGVCFs 对块粒度与 `PL=0,0,0` 块的敏感度）。
- **修复尚未开始**：本报告不含任何代码改动。

## 复现脚本

`fastgatk-native/scripts/verify_gvcf_stream_overlapping_indels_gatk_oracle.py`
（第 1 轮 Track A 产出；当前状态为 **diagnostic**，即报告偏离而不判 pass/fail）。
尚未注册进 CTest；注册需在第 2 轮修复后把断言改为「流式 == 非流式 == GATK」。
