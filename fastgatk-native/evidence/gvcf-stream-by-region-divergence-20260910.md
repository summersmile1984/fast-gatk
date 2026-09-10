# `--stream-by-region` 与 GATK / 非流式输出不一致（2026-09-10）

状态：**已复现，P1（流式路径正确性问题）**。发现于第 1 轮 Track A 审计，
第 2 轮完成定位（**修复尚未开始，本报告不含代码改动**）。

## 结论摘要

`fastgatk-hc-call --stream-by-region N`（当 tile 小于区间长度时）的输出与
**同一条命令去掉 `--stream-by-region`** 的输出不一致，也与 pinned GATK 4.6.2.0 不一致。
第 2 轮已验证出三类不一致：

1. **参考块粒度改变（结构性）**：非流式路径与 GATK 都把 20:10020230–10020428 拆成
   29 个小 `<NON_REF>` 块；流式路径把它合并成 **1 个** 块（`POS=10020230 END=10020428`）。
   在 tile=100/500 时，流式输出比非流式输出**少 28 条记录**。
2. **注解读数错误（数值性）**：在 20:10020680 这一行，流式路径给出的
   `RAW_MQandDP`（`28800,8 → 97200,27`）与 `SB`（`0,0,3,3 → 0,0,0,0`）与非流式、与 GATK 都不同。
3. **凭空产生 phasing（语义性，第 2 轮新发现）**：在 20:10020228/10020229，非流式与 GATK 输出
   未定相的 `GT=0/1`，流式输出为已定相 `0|1`/`1|0` 并新增 `PGT/PID/PS`。

判定基准：非流式路径与该区间上的 GATK **全部 48 条数据行逐字段相同**，所以非流式是对的、
流式是偏离方。第 2 轮另已确认 **只有在 tile 覆盖整个输入区间时才逐字节相同**，
且加大 `--assembly-region-padding` **不是修复**（见「已验证」一节）。

## 重要更正（第 3 轮）：D2 不是流式缺陷，而是 `-L` 窗口依赖缺陷

本报告先前的定位框架是「非流式正确、流式偏离」。**该框架对 D2 是错的**，由独立审计发现、
并经我复核确认：

```bash
# 完全不加 --stream-by-region
fastgatk-hc-call -R fixtures/chr20/ref20mnp.fasta -I fixtures/chr20/mnp.bam \
  -L 20:10020381-10020710 --emit-ref-confidence GVCF --max-mnp-distance 1 --threads 1 \
  --add-output-vcf-command-line false -O n.g.vcf
```

在该窗口下，**没有任何流式**，native 仍然给出错误的
`RAW_MQandDP=97200,27` / `SB=0,0,0,0`。关键对照（同一夹具，两条 `-L`）：

| `-L` 窗口 | GATK | native | GATK vs native 数据行差异 |
| --- | --- | --- | --- |
| `20:10019901-10020710` | `28800,8` / `0,0,3,3`（48 条） | `28800,8` / `0,0,3,3`（48 条） | **0 行（逐字段相同）** |
| `20:10020381-10020710` | `28800,8` / `0,0,3,3`（13 条） | `97200,27` / `0,0,0,0`（13 条） | **2 行** |

即：**GATK 的该注释与 `-L` 窗口起点无关，native 的却有关**。因此

- D2 的真实性质是 **`-L` 窗口依赖**（interval-dependent annotation evidence），
  而非流式路径特有的错误；流式只是通过选择 tile 窗口**触发**了它。
- 这是一个**比流式问题更广的 GATK parity 缺陷**：普通非流式 HaplotypeCaller 只要用不同的
  `-L` 窗口（例如 scatter 分区）就会复现。本报告原标题中的「与 GATK 不一致」对 D2 而言，
  根因不在 streaming。
- 我先前「非流式是正确的一方」这一判断**只对 `20:10019901-10020710` 这一个窗口成立**，
  不能外推。

D1（凭空 phasing）与 D3（参考块粒度）的流式相关性见
`fastgatk-native/evidence/2026-09-10-parallel-audit/track-d-streaming-rootcause.md`：
该审计把 D1 钉到 `gvcf()` 的 phasing 块（`hc_call.cpp:4085-4195`），
把 D3 钉到 tile 内部的逐 locus RCM 取值（tile=500 在 10020230 算出
`depth=8 informative=0 gq=0 pl=0,0,0`，而非流式为 `depth=12 informative=12 gq=36 pl=0,36,494`）。

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

### 更正：先前「证伪拼接器」的推理有误

先前结论建立在**错误前提**上：我假设各 tile 产出的块与非流式的块相同，于是用非流式的 29 个块
去预测 `merge()`（`hc_call.cpp:5864`，PL 取逐元素最小值）的结果，得到
`PL=0,26,494 / GQ=26 / MIN_DP=12`，与实测的 `PL=0,0,0 / GQ=0 / MIN_DP=21` 不符，
据此判断「拼接器不是成因」。

该前提不成立：**各 tile 自己产出的块就已经是退化的（`PL=0,0,0`）**，两个退化块经 `merge()`
自然得到 `PL=0,0,0 / GQ=0`，与实测一致。而且从坐标上可判定拼接器**必然参与**：记录按 core
坐标输出，tile=500 的 core 边界在 10020401，而该记录为 `10020230–10020428`，跨越边界，
单个 tile 无法独立产出它。

正确表述：**拼接器是把两个已退化的块接起来的机制；真正的缺陷是 tile 在 core 边界附近算出的
参考置信度本身是错的**（DP 12–13 而非 31–32、`PL=0,0,0`）。

### 已验证：偏离集中在 tile 边界附近，且随 halo 变化（padding 不是修复）

同一夹具上做 tile × `--assembly-region-padding` 扫描，统计与非流式输出的差异行数：

| `--stream-by-region` | pad=100 | pad=200 | pad=300 |
| ---: | ---: | ---: | ---: |
| 300 | 1 行 | 1 行 | 1 行 |
| 405 | 43 行 | 43 行 | 3 行 |
| 500 | 43 行（仅 20 条记录） | 3 行 | 3 行 |
| 810（= 区间长度，单 tile） | **0 行（逐字节相同）** | 1 行 | 1 行 |

两个结论：

1. **默认 pad=100 下，只有单 tile（tile ≥ 区间长度）才与非流式逐字节相同**；任何真实分块都会偏离。
2. **加大 padding 不是修复**：它能把 tile=500 的块结构修回来（43 → 3 行），却把原本逐字节相同的
   单 tile 情形弄坏（0 → 1 行）。偏离只是被搬到了别处。

### 已验证：多 tile 下有两类不可消除的残留偏离（tile=500, pad=300）

**(D1) 流式路径凭空产生 phasing —— 本次新发现。** 非流式（= GATK）在 10020228/10020229 输出未定相：

```text
10020228  ...  0/1:4,8,0:12:99:302,0,129,314,153,467:3,1,5,3
```

流式输出为已定相，并新增 PGT/PID/PS：

```text
10020228  ...  0|1:4,8,0:12:99:0|1:10020228_G_A:302,0,129,314,153,467:10020228:3,1,5,3
```

`GT` 由 `0/1` 变为 `0|1`/`1|0` 并附带 PGT/PID/PS —— 这是**输出语义**的改变，不只是数值偏差。

**(D2) 注释 read 作用域错误。** 同前：10020680 处
`RAW_MQandDP` 由 `28800,8` 变为 `97200,27`，`SB` 由 `0,0,3,3` 变为 `0,0,0,0`。

进一步测量表明 **D2 是「分块相关」的，而非某个固定 tile 边界的局部效应**。
同一 locus 在不同分块下取值如下（默认 pad=100）：

| `--stream-by-region` | 10020680 的 `RAW_MQandDP` | `SB` | 与非流式一致 |
| ---: | --- | --- | --- |
| 非流式 | `28800,8` | `0,0,3,3` | 基准 |
| 200 | `97200,27` | `0,0,0,0` | 否 |
| 300 | `97200,27` | `0,0,0,0` | 否 |
| 405 | `28800,8` | `0,0,3,3` | 是 |
| 500 | `28800,8` | `0,0,3,3` | 是 |
| 600 | `97200,27` | `0,0,0,0` | 否 |
| 700 | `28800,8` | `0,0,3,3` | 是 |
| 810 | `28800,8` | `0,0,3,3` | 是 |

并且 `--assembly-region-padding` 也会翻转这个结果（tile=500 在 pad=100 正确、pad=300 错误）。
即：**流式路径的注释不满足「分块不变性」，而非流式路径按构造满足。**

已定位到的可疑代码位置（尚未证实，供修复用）：`run_region_streaming` 的每 tile 渲染调用

```cpp
auto result = fastgatk::calling::run(decoded.reads, references, tile_options);   // hc_call.cpp:6666
...
? gvcf(metadata_reader, result, ..., &decoded.reads,                              // hc_call.cpp:6674-6681
       tile_options.informative_read_overlap_margin, options.max_alternate_alleles)
```

`decoded.reads` 是 **tile 的 halo 读取批**，`tile_options` 携带该 tile 的
`interval_start/end` 与 `informative_read_overlap_margin`。注释证据集看起来由此二者共同决定，
因此随分块变化。注意不能简单归因为「batch 过大」：tile=810 的 batch 更大，却给出正确的 `28800,8`。

### 仍未完成：定位到具体代码

候选仍是参考块的两个构造器（精确的 `build_reference_blocks` vs 近似的
`build_profile_local_reference_blocks`），以及 `run_region_streaming` 的
`input_halo_intervals`（`hc_call.cpp:6463`）如何决定 tile 的读取窗口。

`FASTGATK_DEBUG_GVCF_EMISSION=1` 的输出显示 `symbolic-pl=unavailable` 在**两种模式下都出现**，
因此它**不是**判别依据（先前曾把这条列为优先线索，现更正）。

### 与状态文档旧记录的冲突

`IMPLEMENTATION_STATUS.md` 曾记录「流式 gVCF 输出路径已将实际 read batch 传给注释重算，
而非使用空 read 指针」，并据此视为已修正。本次结果表明：该改动**可能传入了非空但范围错误的
read batch（分区级而非 locus 级）**，因此"不再用空指针"并不等于"数值正确"。

## 已证明 / 未证明

**已证明**

- 非流式 native 与该区间上的 GATK 4.6.2.0 **全部 48 条数据行逐字段相同**（`diff` 无输出），
  包括 10020231–10020428 的 29 块碎片结构。这是本报告的基准：非流式路径是正确的一方。
- 在 `fixtures/chr20/mnp.bam` + `20:10019901-10020710` 上，`--stream-by-region N`
  （N < 区间长度）与非流式输出不一致，且**默认参数下不存在任何逐字节相同的多 tile 配置**。
- 偏离随 tile 与 `--assembly-region-padding` 变化：加大 padding 能修回块结构，却会弄坏
  单 tile 情形（见扫描表）——**padding 不是修复**。
- **两个可分别复现的残留缺陷**（tile=500, pad=300）：
  D1 在 10020228/10020229 凭空产生 phasing（`GT 0/1 → 0|1` 并新增 PGT/PID/PS）；
  D2 在 10020680 注释 read 作用域错误（`RAW_MQandDP`、`SB`）。
- OpenMP 与 Serial 的流式偏离完全相同（Host 路径问题，非 kernel）。
- 覆盖范围本身未丢失：tile=500 时 28 条被合并进 `10020230–10020428` 单块，
  10 条变异记录逐一保留。
- **更正**：先前「拼接器 `can_merge` 已被证伪」的结论作废（推理前提有误，见「更正」一节）。
  拼接器必然参与跨 core 边界的块合并；真正的缺陷在 tile 侧的参考置信度取值。

**未证明**

- 未定位到具体代码行；`build_profile_local_reference_blocks` 与
  `build_reference_blocks` 的取舍、以及 `input_halo_intervals` 的窗口口径都只是候选，
  未证实缺陷落在哪一处。
- **D1（凭空 phasing）与 D2（注释作用域）是否同源，尚未判定。**
- `RegionGvcfStitcher::can_merge` 只比较 `pl.size()` 而不比较 PL 值：本例中两个输入块
  都是 `PL=0,0,0`，故未暴露；在其他夹具上仍可能造成不正确的块合并，需单独验证。
- 未在该区间之外确定偏离边界（哪些夹具会触发、是否所有多重叠 indel 窗口都触发）。
- 未验证 `--stream-by-contig`、以及与 `--max-alternate-alleles`、read-overlap margin
  等选项的组合。
- 未验证下游影响（如 GenotypeGVCFs 对块粒度、`PL=0,0,0` 块、以及意外 phasing 的敏感度）。
- **修复尚未开始**：本报告不含任何代码改动。

## 复现脚本

`fastgatk-native/scripts/verify_gvcf_stream_overlapping_indels_gatk_oracle.py`
（第 1 轮 Track A 产出；当前状态为 **diagnostic**，即报告偏离而不判 pass/fail）。
尚未注册进 CTest；注册需在第 2 轮修复后把断言改为「流式 == 非流式 == GATK」。
