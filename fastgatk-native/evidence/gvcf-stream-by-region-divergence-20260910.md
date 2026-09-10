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

## 根因假设（未验证，供第 2 轮定位）

1. **注释 read batch 取错范围**：流式路径把 tile 范围的 read batch（27 条）交给了
   locus 级注释重算，而非与该 locus 重叠的 reads（8 条），于是 `RAW_MQandDP` 反映的是
   tile 统计；`SB` 退化为全 0，可能同源于所传 reads 缺少该 locus 的链/碱基计数。
   这与 `IMPLEMENTATION_STATUS.md` 记录的「流式 gVCF 输出路径已将实际 read batch 传给
   注释重算，而非使用空 read 指针」改动直接相关——该改动可能传入了**非空但范围错误**的 batch。
2. **参考块边界由 tile 作用域决定**：流式写入器在 tile 内把相邻 `<NON_REF>` 块合并，
   而非流式写入器按 locus 级 active 状态切块；tile 边界与 active 状态不一致时块边界丢失。

两者都位于 Host 路径（两个后端输出完全相同，排除 kernel 差异）。

## 已证明 / 未证明

**已证明**

- 在 `fixtures/chr20/mnp.bam` + `20:10019901-10020710` 上，`--stream-by-region N`
  （N < 区间长度）与非流式输出不一致：记录数最多少 28 条；tile=300 时字段
  `RAW_MQandDP`/`SB` 不同。
- 非流式 native 与该区间上的 GATK 4.6.2.0 在记录数与 20:10020680 行上一致。
- OpenMP 与 Serial 的流式偏离完全相同（Host 路径问题）。
- 覆盖范围本身未丢失：tile=500 时 28 条被合并进 `10020230–10020428` 单块，
  10 条变异记录逐一保留。

**未证明**

- 未在该区间之外确定偏离边界（哪些夹具会触发、是否所有多重叠 indel 窗口都触发）。
- 未定位到具体代码行；上述根因仅为假设。
- 未验证 `--stream-by-contig`、`--stream-by-region` 与 `--max-alternate-alleles`、
  read-overlap margin 等选项组合。
- 未验证下游影响（如 GenotypeGVCFs 对块粒度的敏感度）。

## 复现脚本

`fastgatk-native/scripts/verify_gvcf_stream_overlapping_indels_gatk_oracle.py`
（第 1 轮 Track A 产出；当前状态为 **diagnostic**，即报告偏离而不判 pass/fail）。
尚未注册进 CTest；注册需在第 2 轮修复后把断言改为「流式 == 非流式 == GATK」。
