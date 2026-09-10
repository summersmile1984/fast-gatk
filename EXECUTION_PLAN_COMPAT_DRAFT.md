# fast-gatk 兼容性与验证执行方案（草案）

> 本文是把运行时规范、模块实施蓝图和当前 native 原型收敛成可执行的工程验收方案。
> 它定义“可以直接替换 GATK”的证据链，而不是把能解析 `-I/-R/-L/-O` 当成兼容完成。
> 目标运行环境是现有的 `GATK + SLURM + Nextflow`；native 侧统一使用 C++ Host、HTSlib
> I/O、Kokkos Views/ExecutionSpace/parallel primitives。未通过本方案相应门槛的工具，
> registry 状态必须保持 `fallback` 或 `adapter`。

相关基线：

- [FAST_GATK_RUNTIME_SPEC.md](</home/turing-agents/Documents/fast-gatk/FAST_GATK_RUNTIME_SPEC.md>)：资源、I/O、Strict/Fast、队列和 Tool contract。
- [MODULE_IMPLEMENTATION_PLANS.md](</home/turing-agents/Documents/fast-gatk/MODULE_IMPLEMENTATION_PLANS.md)：各计算模块的 Host/kernel 边界、依赖和阶段。
- [GATK_MODULE_INVENTORY.md](</home/turing-agents/Documents/fast-gatk/GATK_MODULE_INVENTORY.md>)：完整工具/源码区域清单。
- [GATK_KOKKOS_RESEARCH.md](</home/turing-agents/Documents/fast-gatk/GATK_KOKKOS_RESEARCH.md)：替换策略、集群契约和 GATK oracle 来源。

## 1. 兼容性定义和不可妥协项

### 1.1 三个兼容等级

每个 tool、backend 和版本在 registry 中明确登记下列状态，不使用模糊的“兼容”：

| 等级 | 含义 | 可否替换原 GATK process |
|---|---|---|
| `byte-identical` | 在固定 toolchain、输入、参考、interval、seed、压缩参数和 Strict 模式下，主输出及规定 sidecar 字节一致 | 可以，但仍须通过 pipeline/资源门槛 |
| `bio-identical` | VCF/GVCF 经过规范化后，记录、等位基因、GT/过滤和规定数值容差与 GATK 一致；允许 header、压缩 block 或浮点低位不同 | 仅限在 tool registry 明确批准的模式 |
| `contract-compatible` | CLI、退出码、文件格式、header/dictionary、索引、sidecar、排序、scatter/gather 和恢复契约一致；算法结果尚未达到 GATK 等价 | 只能作为试运行或显式 native prototype |

`fast` 模式不是隐含的兼容模式：它可以改变浮点归约顺序和 backend，但必须在 manifest
写出 `determinism=fast`、误差界限和 oracle 版本。Strict 模式固定归约树、排序 tie-break、
随机种子、数学函数和 writer 参数。

### 1.2 直接替换的边界

替换发生在进程边界：Nextflow 的 process 名称、输入 channel、输出文件、cache key、
`-resume` 和 scatter/gather 不改变。native 进程不能偷偷把多个原 GATK process 融合，
不能在作业内自行 `sbatch`，不能启动第二套 MPI/Spark 调度。显式 fused profile 只能作为
后续新模式，不能改变 drop-in profile 的语义。

launcher 负责：

1. 识别 launcher 参数、`@args` 文件、`--` 分隔符和 tool 名；保留原始 argv。
2. 查 `tool_registry.json`：tool、版本、参数子集、输入输出类型、sidecar、backend、determinism 和 golden corpus。
3. 按资源、输入类型和 backend 选择 native；不满足契约时显式 fallback 或 fail-closed。
4. 记录选择结果、fallback 原因、版本、输入/reference/index digest 和资源预算。

native tool 不应再次调用 launcher。JVM 专用选项在 native 路径被记录并忽略，在 fallback
路径原样传递；未知参数不能静默丢弃。CUDA/HIP/SYCL 初始化失败是 `BACKEND_UNAVAILABLE`
或 CPU native fallback，不得伪装成 `BAD_INPUT`。

## 2. 当前代码审计结论

### 2.1 已有可复用实现和证据

| 区域 | 当前实现 | 证据/命令 | 结论 |
|---|---|---|---|
| Kokkos 公共边界 | `fastgatk-core` 的 `HostBatch`、`DeviceBatch<ExecSpace>`、`KernelPlan<ExecSpace>` | PairHMM 和四类 demo 已链接同一 target | 架构边界成立；尚未覆盖完整生产 schema/版本化 batch |
| PairHMM | Kokkos SIMD scalar/AVX2/AVX-512 | `python3 pairhmm-demo/scripts/verify_kokkos.py --records=512 --workload=matrix8 --variants=scalar,avx2,avx512`，4096 values `bit_different=0` | CPU Strict 样例通过；尚未覆盖 GATK 全部 PairHMM vectors、flow、GPU |
| SW/BQSR/assembly/region | 规则 DP、整数表、k-mer degree、region stream | `python3 kokkos-modules-demo/scripts/verify_modules.py --records=512 --threads=4`，12 组合通过 | 仅 primitive 验证；SW CIGAR/traceback、真实 BQSR/graph 尚未完成 |
| HTSlib | 真实 BAM/CRAM reader，CRAM 可指定 reference | `python3 fastgatk-native/scripts/verify_native.py` | BAM/CRAM 读取和 batch boundary 通过；VCF/BCF/index writer 仍需独立闭环 |
| CountReads/FlagStat | C++ Host + HTSlib bounded stream，Kokkos RangePolicy（FlagStat 为 12-counter atomic）和 OutputManifest | `python3 fastgatk-native/scripts/verify_count_reads.py`、`python3 fastgatk-native/scripts/verify_flag_stat.py`、`python3 fastgatk-native/dispatcher/verify_dispatcher.py` | pinned chr17 fixture、重复 `-I`、interval/filter/manifest、OpenMP/Serial 与 GATK 4.6.2.0 oracle 通过；registry/binary status 为 `contract-compatible`，remote/cloud 和未覆盖 GATK filter 仍显式 fallback |
| HaplotypeCaller 原型 | `fastgatk-hc-call`：pileup → 单 SNP candidate → 简化 likelihood/genotype | BAM 493 reads、CRAM 8 reads；跨 batch signature 稳定 | 仅 `contract-compatible` prototype，不是 HC 等价 |
| CLI | `fastgatk-native/dispatcher/fastgatk`、`gatk_compat.sh`、registry、`@args`、`-I/-R/-L/-O`、`--help/version/list/dry-run/fallback` | `python3 fastgatk-native/dispatcher/verify_dispatcher.py`、`bash fastgatk-native/scripts/verify_pipeline_local.sh` | P0 dispatcher/registry 基础已通过；仍未覆盖全量 Barclay 参数、完整 GATK launcher 行为 |
| VCF/manifest | VCF/VCF.GZ、manifest、telemetry | `verify_native.py` 和 local pipeline verifier | 文件边界可用；尚未生成 VCF index，manifest 尚未做完整排序/interval/sidecar 验证 |
| GATK oracle | `verify_gatk_oracle.py` 同一 chr17 区域调用 GATK 4.6.2.0/native | sentinel 共享，call-set intersection=1；manifest 明示 `bit_identical=false` | oracle 通路成立，但算法差异很大，不能宣称替换 |
| Nextflow/SLURM | DSL2 smoke、SLURM 本地模拟 wrapper | local verifier；真实 Nextflow 仅在安装时执行 | process 契约样例通过；真实 SLURM scatter/gather/resume/retry 尚未测 |

### 2.2 当前必须保持的诚实边界

当前 HC/Mutect2 调用链为最小单 SNP/局部 somatic prototype，缺少 CIGAR-aware/indel assembly、真实
haplotype PairHMM、GATK genotype priors、reference confidence/GVCF、VCF index、
完整 header/PG/annotation、Mutect2 sidecar、BQSR report 和 GenomicsDB。P0
dispatcher/registry 已能明确路由 native、打印 dry-run 计划，并对未知参数返回
结构化错误或在显式 `--fallback` 时转交配置的 Java 命令；它还不是全量 Barclay
launcher。CountReads/FlagStat 已满足本草案定义的 CLI、退出码、文本输出、输入
遍历、过滤、确定性和 manifest contract gate，可标成 `contract-compatible`；该等级
不等于所有 GATK 工具均已实现，也不代表 HC/Mutect2 的算法结果等价。

## 3. 目标代码与产物布局

生产仓库应固定为以下依赖方向：

```text
compat/       gatk launcher, @args, registry, fallback
runtime/      cgroup/SLURM probe, controller, queues, telemetry
traversal/    intervals, shard/core+halo, deterministic order
io/           HTSlib, BGZF/CRAM/VCF/BCF, reference, spill, indexes
model/        flat SoA schemas, arena, versioned batch/chunk
kernels/      PairHMM, SW, k-mer, pileup, BQSR, genotype, dense math
tools/        HC, Mutect2, BQSR, GenotypeGVCFs, CNV, VQSR adapters
formats/      OutputManifest, schema, error codes, VCF/BAM metrics
tests/        kernel, format, tool, oracle, pipeline, resource, fuzz
```

一个发布 artifact 必须同时包含：

- `fastgatk`/`gatk` 兼容 launcher；
- native binaries 及支持的 Kokkos backend；
- `tool_registry.json`、schema version、build/compiler/Kokkos/HTSlib provenance；
- golden corpus manifest（输入/reference/index digest、预期输出 digest/规范化 digest）；
- `compatibility-report.json` 和版本化 telemetry schema；
- Java GATK fallback 入口及 fallback reason code。

## 4. CLI 和 launcher 替换执行方案

### 4.1 解析顺序

实现 `fastgatk` dispatcher，保持以下行为：

```text
launcher options → @args expansion → tool name → tool options → -- separator
```

必须覆盖：

- `--help`、`--version`、`--list`、`--dry-run`、`--java-options`、`--gatk-config-file`；
- `@file` 递归深度上限、空白/引号/转义规则和 argv 原顺序；
- GATK 常用短/长别名：`-I/--input`、`-R/--reference`、`-L/--intervals`、`-O/--output`；
- tool-specific 参数不由 launcher 重新解释；原始 argv 必须可供 fallback；
- 参数缺失、重复冲突、未知参数和输入不存在时的 exit code/错误分类。

`gatk_compat.sh` 是当前最小入口，只识别部分 launcher 参数。它应作为 dispatcher
合同测试的前向兼容样例，不能直接当作生产入口。

### 4.2 Registry 条目

每个 tool entry 至少包含：

```json
{
  "tool": "HaplotypeCaller",
  "gatk_version": "4.6.2.0",
  "native_binary": "fastgatk-hc-call",
  "status": "fallback|adapter|native",
  "supported_parameters": [],
  "required_outputs": ["primary", "index", "stats", "sidecars"],
  "backends": ["openmp", "cuda"],
  "determinism": ["strict", "fast"],
  "golden_corpus": "corpus/hc-v1.json",
  "resource_model": "hc-region-v1",
  "fallback_policy": "explicit|fail_closed"
}
```

entry 的状态晋级必须由第 10 节的 gate 驱动，不能由 binary 存在或 CLI 解析成功驱动。

### 4.3 CLI 验收

```bash
set -euo pipefail
root=/home/turing-agents/Documents/fast-gatk

FASTGATK_HC_BINARY="$root/fastgatk-native/build/fastgatk-hc-call" \
  bash "$root/fastgatk-native/scripts/gatk_compat.sh" --java-options '-Xmx1g' \
  HaplotypeCaller -I "$root/gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam" \
  -L 17:69000-69100 -O /tmp/fastgatk-calls.vcf

"$root/fastgatk-native/build/fastgatk-hc-call" --help
bash -n "$root/fastgatk-native/scripts/gatk_compat.sh"
```

生产 CI 还必须加入 `@args`、等号形式、`--`、坏参数、未知 tool、fallback 和
`--dry-run` 的 golden stderr/exit-code 检查。

## 5. 输入、输出和格式契约

### 5.1 输入

- BAM/CRAM/SAM 由 HTSlib 读取；CRAM 必须显式验证 reference、FAI、sequence dictionary。
- VCF/BCF、tabix/CSI、FASTA/FAI、dict 的 header/contig 顺序和 digest 写入 run metadata。
- remote URI、range request、cache/spill 由 Host 管理；kernel 不接触文件或网络。
- 每个 decoded chunk 携带 `job_id`、`sample_id`、`shard_id`、`core_interval`、`halo_interval`、`input_offset`、`output_ordinal`、`schema_version`、CRC。

### 5.2 输出

native writer 关闭前必须完成：

1. primary output 存在、可读、非 partial，VCF/BCF/BAM/CRAM 的 header 和 contig dictionary 完整。
2. coordinate/sample/allele/order 合法；core interval 无重复、无漏写，halo 不输出。
3. 索引类型正确（BAI/CSI/CRAI/TBI），并能被 HTSlib 重新打开、按 region 查询。
4. 规定的 `.stats`、F1R2、metrics、recalibration report、gVCF/sidecar 全部存在。
5. OutputManifest 记录 path、kind、size、SHA-256、schema、complete、interval set、index、producer、fallback、determinism 和 telemetry。
6. 任一验证失败时删除或标记 partial output，退出 `OUTPUT_CONTRACT_FAILURE`，不能返回成功。

默认输出不引入私有格式；spill/run 文件应有版本化 header、长度、CRC、排序键、shard
ID 和输入 digest，成功后按 manifest 清理，失败后保留最后一个完整 run 供恢复。

### 5.3 格式验收命令

```bash
root=/home/turing-agents/Documents/fast-gatk
python3 "$root/fastgatk-native/scripts/verify_native.py"
python3 "$root/fastgatk-native/scripts/verify_count_reads.py"
python3 "$root/fastgatk-native/scripts/verify_flag_stat.py"
"$root/fastgatk-native/build/fastgatk-count-reads" --help >/dev/null
"$root/fastgatk-native/build/fastgatk-flag-stat" --help >/dev/null
bash "$root/fastgatk-native/scripts/verify_pipeline_local.sh"

# 在安装 samtools/bcftools 的 CI image 中增加：
samtools quickcheck -v calls.bam
samtools view -H calls.bam
samtools index -c calls.bam
bcftools view -h calls.vcf.gz
bcftools index -t calls.vcf.gz
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' calls.vcf.gz
```

CRAM、VCF.GZ 当前已做读写 smoke；真实生产 gate 还要增加 CRAI、CSI/TBI、BCF、gVCF
和 header/PG/index read-back。不能用“VCF 能打开”替代这些检查。

## 6. Determinism、bit-identical 和 numerical policy

### 6.1 数据层确定性

Strict 模式固定：

- 输入 read 的 coordinate/read-name/order、downsampling seed、reference/contig/allele/sample order；
- length bucket、pair permutation、graph key、candidate stable ID、scatter output ordinal；
- per-thread/per-team 局部表的合并顺序和 reduction tree；
- NaN/underflow/zero、Phred/error table、log10/exp 实现、`-ffp-contract=off` 等编译契约；
- BGZF/VCF writer 的 compression level、header/PG 行顺序和 index 生成参数。

Fast 模式允许 device/BLAS 的浮点顺序差异，但必须输出 `tolerance_profile`，并通过
规范化 VCF、数值误差、call-set 和 truth-set gate。Strict 失败不能被 Fast 结果掩盖。

### 6.2 验收分层

| 层 | Strict gate | Fast gate |
|---|---|---|
| PairHMM | matrix/vector raw bits；GATK table/log10 固定 | likelihood abs/rel tolerance、GT/GQ 和 call-set |
| SW | score、CIGAR、start/end、overhang、tie-break exact | score exact 或限定差异；CIGAR 必须相同/可规范化 |
| BQSR | 整数 table/report、read order exact | table exact；仅并行 telemetry 可不同 |
| assembly | stable graph/path/haplotype order | canonical graph/path 和 truth-set 一致 |
| genotyping | PL/AD/DP/GQ 按工具指定 raw/tolerance | GT/filter/call-set/truth gate |
| VCF/GVCF | header/record/index/sidecar byte 或 canonical digest | canonical records、block、index 可读 |

当前证据只有 PairHMM 4096 matrix 的 scalar/AVX2/AVX-512 raw-bit 一致；
`fastgatk-hc-call` 的整数 `signature` 只证明跨 batch/thread 的 prototype 稳定，
不证明浮点 bit identity 或 GATK biological equivalence。

### 6.3 Determinism 验收命令

```bash
root=/home/turing-agents/Documents/fast-gatk
python3 "$root/pairhmm-demo/scripts/verify_kokkos.py" \
  --records=512 --workload=matrix8 --variants=scalar,avx2,avx512
python3 "$root/kokkos-modules-demo/scripts/verify_modules.py" \
  --records=512 --threads=4
python3 "$root/fastgatk-native/scripts/verify_native.py"

# 生产 gate：同一输入重复、换 batch、换线程和重新启动；比较 manifest 的
# canonical_output_digest、record digest、signature、排序和 telemetry schema。
for seed in 1 2 3; do
  FASTGATK_DETERMINISM=strict FASTGATK_SEED="$seed" \
    "$root/fastgatk-native/build/fastgatk-hc-call" ...
done
```

严格测试必须至少覆盖空输入、单 read、同坐标并列 read、N/低质量碱基、最大 read/
haplotype 长度、duplicate、unmapped、跨 contig、边界 interval 和批次恰好切分点。

## 7. Golden corpus 和 GATK oracle

### 7.1 Corpus 分层

建立版本化 `corpus/manifest.json`，每项记录：GATK 版本/JDK、native commit、输入和
reference/index SHA-256、interval、tool args、seed、模式、预期 output digest、
canonical record digest、metrics/sidecar digest、资源配置和允许 tolerance。

四层 corpus：

1. **kernel vectors**：GATK `PairHMMUnitTest`、`VectorPairHMMUnitTest`、SW random/edge vectors、quality/error table；数值 raw/tolerance。
2. **format vectors**：SAM/BAM/CRAM、VCF/BCF、CSI/BAI/CRAI/TBI、header/PG、压缩和坏输入；由 HTSlib/samtools/bcftools read-back。
3. **tool vectors**：GATK `HaplotypeCallerIntegrationTest`、`Mutect2IntegrationTest`、BQSR/GenotypeGVCFs/GenomicsDB/CNV/VQSR representative fixtures；比较完整 primary+sidecars。
4. **production/truth vectors**：NA12878/NA24385 GIAB 小区域、肿瘤/正常样本、真实 scatter shards；用 hap.py/GA4GH/NIST truth/confident regions 评估，不只做文本 diff。

### 7.2 Oracle 对比程序

`verify_gatk_oracle.py` 目前只验证一个 chr17 sentinel 和 call-set overlap，并明确
`bit_identical=false`；它是结构性 smoke，不能成为 release gate。下一版 oracle runner
应为每个 tool 生成：

```text
GATK output → raw file digest
            → header-normalized digest
            → coordinate/allele-normalized records
            → numeric fields (PL/GQ/QUAL) tolerance report
            → sidecar/metrics digest
            → truth-set metrics
```

验收规则：

- 缺少 GATK/JDK/reference/index 时 gate `SKIP` 只能用于非 release PR，release 必须 `FAIL`；
- native 多出/少出记录必须列出 first-difference 和 interval；
- manifest 宣称 `bit_identical=true` 时必须有 raw digest 证据；
- 只共享 sentinel 不能标记 native；
- oracle 版本升级时建立新 corpus version，不覆盖旧基线。

### 7.3 GATK oracle 命令

```bash
root=/home/turing-agents/Documents/fast-gatk
FASTGATK_REQUIRE_GATK_ORACLE=1 \
  python3 "$root/fastgatk-native/scripts/verify_gatk_oracle.py"
```

未来把同一 runner 扩展为：

```bash
./ci/run_oracle.sh PairHMM --corpus corpus/pairhmm-v1.json --mode strict
./ci/run_oracle.sh HaplotypeCaller --corpus corpus/hc-v1.json --mode strict
./ci/run_oracle.sh Mutect2 --corpus corpus/mutect2-v1.json --mode strict
./ci/run_oracle.sh BaseRecalibrator --corpus corpus/bqsr-v1.json --mode strict
```

## 8. 所有模块的执行顺序、状态和 gate

| 模块/工具 | 首版实现策略 | 依赖 | native 晋级 gate |
|---|---|---|---|
| compat/dispatcher | C++ Host；registry、`@args`、fallback、错误码 | 无 | CLI/exit/fallback/registry golden |
| runtime/traversal | C++ Host；cgroup/SLURM probe、core+halo、byte queues、adaptive batch | compat、model | 资源 cap、边界、spill、scatter 无重复 |
| HTSlib/model/formats | C++ Host；BAM/CRAM/VCF/BCF/index、flat SoA、manifest | runtime | htslib/samtools/bcftools read-back |
| PairHMM regular/flow | Kokkos SIMD CPU；Team/Vector GPU；Host table/workspace | model、plan | GATK/GKL vectors，raw/tolerance、性能 |
| Smith-Waterman | Kokkos score/wavefront；Host traceback/CIGAR | model、PairHMM 前置 | score+CIGAR+overhang+tie-break |
| k-mer/read threading/assembly | Host graph + Kokkos encode/sort/count/score | traversal、SW | graph/path/haplotype order、indel corpus |
| pileup/activity/BAQ/downsampling | Host state machine + Kokkos window histogram/filter | traversal、io | seed/order/window boundary/filter oracle |
| HaplotypeCaller | Host orchestration，复用 assembly/SW/PairHMM/genotype | 上述 P0/P1 | HC VCF/GVCF/index/header/sidecar、GIAB |
| likelihood/genotype/reference confidence | Host allele/block state + Kokkos small-matrix math | PairHMM、HC | PL/AD/DP/GQ、GVCF block/canonical digest |
| Mutect2 | 复用 HC kernels；somatic filter/contamination/orientation Host | HC、likelihood | VCF + stats + F1R2 + somatic GVCF |
| BaseRecalibrator | Host covariates；Kokkos local integer tables/merge | io、runtime | report/table/read-order exact |
| ApplyBQSR | Host streaming transform；可选 Kokkos quality batch | BQSR | BAM/CRAM read-back、quality exact |
| GatherBQSRReports/AnalyzeCovariates | Host/report adapter | BQSR | report schema/metrics |
| GenomicsDBImport | 复用 GenomicsDB/TileDB；Host resource/fd/tmp adapter | formats、runtime | workspace read-back、batch/restart |
| GenotypeGVCFs/CombineGVCFs/ReblockGVCF | Host locus/block/writer + Kokkos genotype math | GenomicsDB、genotype | GVCF block/allele/sample/order |
| MarkDuplicates/SortSam | P2 external-memory Host；早期 Picard fallback | io、spill | metrics、duplicate flags、sort/index、restart |
| VQSR | CPU model control/convergence；Kokkos Kernels VBEM/GMM | dense math、formats | model/tranche/report/tolerance |
| CNV/CollectReadCounts | Kokkos interval overlap/histogram；Host HDF5/TSV | io、dense math | counts/HDF5/SVD/segments/restart |
| Funcotator/annotation | CPU cache/index/string；不强求 GPU | formats、data source | VCF/MAF/SEG field/order/header |
| variant utilities/reference/coverage | HTSlib/CPU；可选 batch filters | formats | format/semantic tests |
| SV/BWA/PathSeq/Spark/Picard 长尾 | fallback-first；按 profile 单独立项 | external deps | tool-specific oracle 后才 native |

执行依赖采用以下顺序：

```text
P0 compat + runtime + HTSlib/model
  → P0 PairHMM + SW + k-mer primitives
  → P0 HC file-boundary/tool smoke
  → P1 activity/pileup + full assembly + likelihood/genotype
  → P1 HC/GVCF + Mutect2 + BQSR
  → P2 GenomicsDB/GenotypeGVCFs + external-memory + VQSR/CNV
  → P3 GPU optimization + annotation/SV/BWA/PathSeq/Spark selective native
```

## 9. Nextflow、SLURM、scatter/gather、resume/retry

### 9.1 资源映射

native 在现有 allocation 内读取：

- CPU：`SLURM_CPUS_PER_TASK`，并校验 process affinity；Kokkos worker 不能超过 allocation；
- memory：cgroup v2 `memory.max/current` 优先，回退 `SLURM_MEM_PER_NODE`，不能使用物理机总内存；
- device：`CUDA_VISIBLE_DEVICES`/SLURM GRES 和 backend free memory；
- scratch：`SLURM_TMPDIR` → `TMPDIR` → 配置，检查可写字节和 inode；
- network/remote：读写延迟 EWMA、cache/spill byte budget。

默认目标是 host 60%、device 60%、scratch 80%；两个窗口超过 80% 则减半，低于 50%
且队列空时增加 10–20%。分配失败缩小 batch 只允许重试一次，随后返回
`RESOURCE_EXHAUSTED` 让 Nextflow/SLURM 处理。

### 9.2 scatter/gather

scatter 以 sample/contig/interval/shard 为单位。每个 shard metadata 固化：

```text
shard_id, core_interval, halo_interval, output_ordinal,
input_digest, reference_digest, schema_version, deterministic_seed
```

规则：

1. reader 可读取 halo，kernel 可计算 halo，但 writer 只提交 core interval。
2. gather 按 `output_ordinal` 和 contig coordinate 排序，不按完成时间合并。
3. gather 前验证 core interval 不重复、不漏写、header/dictionary 一致、所有 sidecar 齐全。
4. empty shard 有显式 manifest，不得让 gather 误判为缺失。

### 9.3 resume/cache

Nextflow process 输出文件名和 sidecar 固定；cache key 至少包含 binary/schema/tool
args、input/reference/index digest、interval、backend、determinism、seed 和资源影响
输出的版本。manifest 的 `complete=true`、primary/index/sidecar digest 全部通过才可
被 `-resume` 复用。partial output、旧 schema 或不同 reference digest 一律 cache miss。

### 9.4 retry/failure recovery

错误分类和动作：

| 错误 | 进程动作 | 是否自动 retry |
|---|---|---|
| `BAD_INPUT` | 立即失败，记录第一个坏文件/record | 否 |
| `RESOURCE_EXHAUSTED` | 保留诊断 manifest，Nextflow 用更大 memory/更小 shard retry | 是，次数有限 |
| `BACKEND_UNAVAILABLE` | CPU native 或 Java fallback，记录原因 | 仅策略允许时 |
| `NUMERICAL_CONTRACT_FAILURE` | Strict 失败；禁止静默降级 Fast | 否，需人工/版本修复 |
| `OUTPUT_CONTRACT_FAILURE` | 清理/隔离 partial，重启 shard | 可 retry 一次 |
| `IO_TRANSIENT` | 有界退避；远端 range/cache 重试 | 是，带上限 |
| `FALLBACK_EXECUTED` | 以 fallback 结果继续，但 manifest/日志显式标记 | 不重复 native |

生产必须测试 kill -9、SIGTERM、节点超时、共享盘不可写、scratch 满、inode 满、
device OOM、fd 上限、远端暂时失败和重复执行。最后一个完整 spill run 可恢复，
但不能依赖进程内内存状态。

### 9.5 工作流验收命令

```bash
root=/home/turing-agents/Documents/fast-gatk
bash "$root/fastgatk-native/scripts/verify_pipeline_local.sh"

# 若环境已安装 Nextflow：
nextflow run "$root/fastgatk-native/workflow/nextflow_smoke.nf" \
  --input "$root/gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam" \
  --region 17:69000-69100 \
  --binary "$root/fastgatk-native/build/fastgatk-hc-call" \
  --outdir /tmp/fastgatk-nextflow-run --threads 2 -resume

# 本地模拟 SLURM；真实集群时设置 FASTGATK_USE_SBATCH=1：
SLURM_CPUS_PER_TASK=2 \
  bash "$root/fastgatk-native/workflow/slurm_smoke.sh" \
  "$root/fastgatk-native/build/fastgatk-hc-call" \
  -I "$root/gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam" \
  -L 17:69000-69100 -O /tmp/fastgatk-slurm.vcf
```

真实 cluster gate 需使用至少两个 interval、两个样本、一次失败 retry、一次 resume、
一次 gather，并检查 output ordinal、duplicate core、cache hit 和 SLURM memory/tmp/GRES。

## 10. CI/CD 和发布门槛

### 10.1 PR 快速门

每个 PR 运行：

```bash
set -euo pipefail
bash -n fastgatk-native/scripts/*.sh pairhmm-demo/scripts/*.sh kokkos-modules-demo/scripts/*.sh
python3 -m py_compile fastgatk-native/scripts/*.py pairhmm-demo/scripts/*.py kokkos-modules-demo/scripts/*.py
bash pairhmm-demo/scripts/build_kokkos_variants.sh
bash kokkos-modules-demo/scripts/build_variants.sh
bash fastgatk-native/scripts/build_native.sh
python3 pairhmm-demo/scripts/verify_kokkos.py --records=512 --workload=matrix8 --variants=scalar,avx2,avx512
python3 kokkos-modules-demo/scripts/verify_modules.py --records=512 --threads=4
python3 fastgatk-native/scripts/verify_native.py
python3 fastgatk-native/scripts/verify_count_reads.py
python3 fastgatk-native/scripts/verify_flag_stat.py
bash fastgatk-native/scripts/verify_pipeline_local.sh
```

### 10.2 Nightly/Release 门

- GATK 4.6.2.0 + pinned JDK oracle，全部 kernel/tool corpus；
- x86 generic/AVX2/AVX-512、ARM/OpenMP、CUDA、HIP、SYCL（可用 runner 才执行，release 缺失 required backend 必须标明）；
- `-fsanitize=address,undefined`、浮点异常/NaN、线程 race、OpenMP/Kokkos tools；
- HTSlib/samtools/bcftools format round-trip、压缩/索引、坏输入 fuzz；
- cgroup memory、scratch/inode/fd、device OOM、spill/resume/retry、SIGTERM/kill；
- Nextflow real DSL2 + SLURM allocation + scatter/gather + retry/resume；
- performance regression：decode、prepare、kernel、collect、compress/write、spill 和 end-to-end 分段记录，超过阈值阻止发布。

每个 artifact 生成 provenance：compiler、flags、Kokkos/HTSlib/zlib/JDK/GATK 版本、CPU
feature、backend、container digest、输入/reference/index digest。不能只保存“pass”。

## 11. 风险和缓解

| 风险 | 影响 | 缓解/退出条件 |
|---|---|---|
| GATK/GKL 与 SIMD/GPU 浮点顺序不同 | bit identity 失败、call 差异 | Strict 固定函数/归约；Fast 明示 tolerance；不能用 tolerance 掩盖 Strict |
| assembly/indel 复杂度被低估 | HC 结果差异大 | Host-first graph/path/CIGAR；先完成 oracle 再移植 kernel |
| header/index/sidecar 不完整 | Nextflow gather/cache 下游失败 | OutputManifest 关 writer 前全量验证；任何缺件 fail closed |
| CIGAR/region/halo 边界错误 | 重复/漏计/错误 variant | core-only writer、canonical interval tests、scatter round-trip |
| device OOM/pinned memory 超预算 | 作业被杀且无法恢复 | 统一 byte budget、bounded queues、spill、一次缩 batch 后返回 RESOURCE_EXHAUSTED |
| 共享盘/远端 I/O 成为瓶颈 | native kernel 加速无端到端收益 | reader/writer CPU pool、local scratch、range prefetch、分段 telemetry |
| OpenMP/Kokkos oversubscription | 性能退化/不稳定 | affinity 校验、线程数不超过 SLURM allocation、CI 固定 cores |
| GenomicsDB/HDF5/Picard/Spark 外部语义 | 重写成本和兼容风险 | adapter/fallback-first，不自造 backend；native 仅在 tool-specific oracle 后启用 |
| fallback 造成结果/日志不透明 | 误称 native、排查困难 | manifest `execution_mode`、`fallback_reason`、原始 argv、版本和退出码 |
| current smoke 被误当 production | 业务风险 | registry status/manifest 明示 prototype；release gate 要求 full corpus/truth/pipeline |
| 只在小 fixture 上 benchmark | 错误架构决策 | 加入真实 cohort、长 read、高 coverage、不同 shard/scratch/network；同时看端到端 |

## 12. 里程碑和完成定义

### M0：兼容骨架（当前到下一阶段）

交付 dispatcher、registry、`@args`、错误码、resource probe、OutputManifest/schema、
HTSlib VCF/BCF/index、golden corpus harness。现有 smoke 只作为回归，不升级为 native。

### M1：HC CPU Strict

完成 CIGAR-aware assembly、SW CIGAR/tie-break、真实 haplotype PairHMM、genotyping、
VCF/GVCF/reference confidence、index 和 GATK HC oracle；完成 Nextflow/SLURM scatter,
gather,resume,retry。满足后才允许 `HaplotypeCaller: native(strict)`。

### M2：Mutect2/BQSR/GenotypeGVCFs

完成 sidecar、F1R2/stats、recalibration report、GenomicsDB adapter、gVCF block、
external-memory 和 restart；每个 tool 独立晋级，不因 HC 通过自动继承。

### M3：多 backend 和性能

CPU generic/AVX2/AVX-512、ARM、CUDA/HIP/SYCL correctness；随后以 profile 驱动
Team/Vector/tile 优化。发布必须同时有 kernel-only、amortized 和 end-to-end 证据。

### 完成定义

只有当某个 tool 同时满足以下条件，才可标记 `native`：

1. CLI/exit/fallback 和 GATK process 可直接替换；
2. 输入、主输出、索引、header、排序、sidecar、manifest 均可回读；
3. Strict/Fast 数值模式的 oracle、determinism、truth/call-set 门槛通过；
4. 单机资源上限、scratch/fd/device OOM、spill、retry、resume 通过；
5. Nextflow + SLURM scatter/gather 的 core/halo 和 cache 契约通过；
6. 至少一个真实 cohort/production-like corpus 和性能回归门通过；
7. fallback、版本、reference/index、backend、tolerance 和未实现参数全部透明记录。

任一条件缺失，状态保持 `adapter`/`fallback`/`contract-compatible`，并在 manifest
中明确不能宣称 bit-identical 或 GATK biological equivalence。
