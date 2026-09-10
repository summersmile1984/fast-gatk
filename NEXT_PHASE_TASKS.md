# Fast-GATK 下一阶段任务清单（基线 54% → 目标 ≥65%）

> 生成日期：2026-09-01（更新：2026-09-03）  
> 基线：`progress_score.json` 全局 **0.5385（54%）**（专家加权，不随 registry 机械变化）；registry 48 工具已定稿（47 contract-compatible、1 adapter、0 prototype）；双后端 CTest 基线 **210/210**（+ chr20 真实数据 contract 后 211 项）。  
> 本文件是下一阶段的可执行清单；总体原则、验收层级、性能门禁仍以 [FAST_GATK_EXECUTION_PLAN.md](/home/turing-agents/Documents/fast-gatk/FAST_GATK_EXECUTION_PLAN.md) §7/§8/§11 为准。

## 执行日志
- **2026-09-03（Round 51）**：chr20 真实参考获取配方固化为可复用脚本 `fastgatk-native/scripts/fetch_chr20_reference.py`（LFS Range 抓取 + UCSC 双锚点坐标校正 + 全长度脚手架/dict/fai 合成，--help/语法校验通过）——A3 数据扩展将可低成本产出新真实窗/参考。
- **2026-09-03（Round 49）**：chr20 真实数据 pin 升级为**完整记录断言（含 REF/ALT）**（4 窗 10 calls；位点集 + 确定性之外，任何等位翻转即 fail）；canonical 实测记录值入库（注意：20228/20229/19967 等密集 SNP 位点等位随配置可翻转——A/G/T 均可能为 read 支持等位，本 pin 锁当前确定性行为作回归防护，非"唯一真理"声明）。OpenMP/Serial 1/1 通过。
- **2026-09-03（Round 48）**：最终 211 项双后端全量基线锁定：OpenMP 211/211、Serial 211/211（含 chr20 真实数据 contract，exit 0 双绿）。
- **2026-09-03（Round 43）**：chr20 真实数据 pin 扩至**四个独立单窗**（新增 10020650-10020710：位点 {10020679,10020680,10020681}，alt 与 read 证据核验一致），OpenMP/Serial contract 1/1 通过（10 calls / 4 窗，位点集合 + 确定性断言）。
- **2026-09-03（Round 42）**：维护收尾——`.all_ctests.txt` 刷新至 211 项（原 200 过时）；`NEXT_PHASE_TASKS.md` 头部基线更新为当前真实状态（48 工具定稿 47/1/0、CTest 210/210 + chr20 contract）。
- **2026-09-03（Round 41）**：`verify_all.sh` serial 门禁正则加入 `fastgatk-hc-chr20-real-contract`（OpenMP 分支本就全量含新测试）——真实数据 pin 纳入最终回归门禁。
- **2026-09-03（Round 40）**：产出 `SESSION_HANDOFF.md` 会话交接文档（已达成/已知限制/资产/建议下一步/纪律），最终门禁复核绿（dispatcher=0、progress pass、hc 子集 3/3）。
- **2026-09-03（Round 39）**：chr20 真实 fixture pin 扩展至**三个独立单窗**（新增 10020390-10020480：位点 {10020429,10020431,10020434}，其 alt 与 reads 全一致 C>G/T>A/A>T），OpenMP/Serial `fastgatk-hc-chr20-real-contract` 1/1 通过（共 7 calls/3 窗，位点集合 + 字节确定性断言）。A3 真实数据覆盖：3 个真实窗口位点锁定（read 证据核验过的确定行为基线）。
- **2026-09-03（Round 38）**：**第二个真实 fixture 固化为 CTest**：新增 `fastgatk-native/scripts/verify_hc_chr20_real_contract.py` + CTest `fastgatk-hc-chr20-real-contract`（存在性条件 add_test），对真实 chr20 mnp.bam 两个**独立单窗**（19901-19970、20200-20500）断言：call 位点集合精确（{10019967,10019969}、{10020228,10020229}）+ 输出字节确定性 + reads≥10。设计说明：等位身份在密集 SNP 簇随窗几何敏感（native 与 GATK 4.6.2.0 均可能翻转），故只 pin 位点集合并记录 read 证据于档案；多 -L/连续窗语义属 A-line 已知限制不纳入。OpenMP/Serial 新 CTest 1/1 通过。扫描 7 个 oracle 无可用的小数值间隙可收紧（Mutect2 TLOD delta 3.5 属 A4 实质差距）。
- **2026-09-03（Round 37）**：canonical 全量双后端基线复核并锁定：OpenMP **210/210**、Serial **210/210**（回退后 exit 0 双双通过）。`PROGRESS_AUDIT_2026-08-31.md` 顶部补"A-line 已知限制"权威记录（缺陷内容/等位错误证据/修复需架构级重构/修复前多 region 连续窗输出以 GATK+read 复核为准）。fixtures/chr20 收编为最小自洽复现集（ref20mnp.fasta+dict+fai、mnp.bam+bai、gatk/native VCF 与 native.json），108MB。执行日志与此前 35 轮证据链一致。
- **2026-09-03（Round 36）**：**A-line 多 region 缺陷修复中止并整体回退（决定性结论）**。深度 2048 + 大 cap 的连续区间实验暴露出比"漏检"更严重的问题：连续单 -L 下**等位身份错误**（20228 真实 ALT=C 却输出 G>A；20434 真实 T 输出 A>G；20438 真实 C 输出 G），而分离 interval/隔离窗均输出正确等位——连续链式多窗口语义（窗拼接/等位合并）尚未正确，继续修复将引入错误 call 回归。**决策**：回退 Round 33-35 全部实验性修复（per-region assemble、窗跨度界、深度分窗、断言改动）至 canonical，恢复原始语义与断言（hc 子集 5/5、双后端重建绿）。全部诊断证据与代码地图归档本文件（R18-36），缺陷作为 **A-line 已知限制**记录：连续多 AssemblyRegion 窗下的 call 完整性与等位正确性需一次专门的架构级重构（逐 region 独立图/装配/基因型 + 窗间等位合并语义），此前 fixtures/chr20 单窗/分离 interval 行为仍与 read 证据及 GATK 一致。
- **2026-09-03（Round 35）**：判别实验——**三分离 interval（A/中 B/远 C）全部 region 正常出 calls**（19967/69、20228/29、20680/81；groups=3、cands 149、hap 144）→ **多 group 下游/partitioned 本身无碍**（此前归因撤回）。残余缺陷精确限定为**连续链式 calling_regions 几何**（full 窗 19901-21411，groups=5，20228/29/20429-38 落在链内部 regions 仍 0/0，候选/单体型/边际化齐备）；疑与连续链分窗接缝/子窗种子/边界 locus 归属相关。修复①②+深度分窗已落地（两/三 interval 均完整出 calls）。子代理 900ba714 带新边界继续定位；回归子集 8/8 绿。
- **2026-09-03（Round 33）**：**双缺陷两处修复落地**。① assemble 候选预算改 **per-calling-region**（regions>1 无条件分区满预算拼接，`calling_pipeline.cpp` run() 7237 起）；② **图窗合并加跨度上限**（= coalesce 的 maximum_span=region+2·halo，阻断跨 gap 链式合并成单种子长窗）。两处均按 GATK AssemblyRegion 语义；同步把 2 个依赖"全局预算抑制必然发生"的 telemetry 断言（`verify_hc_broad_gatk_oracle.py`、`verify_hc_variant_annotations.py` 的 `candidate_softclip_suppressed>=1`）更新为计数存在性断言（附语义注释）。OpenMP/Serial 5/5 子集绿（broad/annotations/kmer-list/streaming/mutect2）。**残余**：full 窗 chr20 中间簇（20228/20229/20429-38，候选与单体型现齐备）仍无 call，仅首尾区 calls——多 group 下游候选→call 映射缺陷（子代理 900ba714 定位中）；修复后再验 full 窗 + 双后端 210/210 + 防复发 CTest。
- **2026-09-03（Round 32）**：大 cap 对照实验（决定性）：full 窗 `--max-candidates 100000` → candidate_sites 1072、calls = 19967(C>G)/19969/20229/20680 共 5——**20228（24/24 reads C）与 20429-38 簇（31-32 reads 全 alt）即使候选海量仍缺失**，而相邻 20229 却出现 → **①全局 64 cap 已排除为唯一成因**；存在第二层缺陷：特定真实位点（20228、20429 簇）在合并长图窗下**无本地单体型/基因型化覆盖**（单种子 DFS 枚举路径未覆盖；reads 支持度极高故非 pruned）。此项已完整归档为 **A-line 开放实现缺陷**：需逐 region 独立图+装配+PairHMM+基因型的重构（GATK AssemblyRegion 语义），修复后以 fixtures/chr20 全窗含 20228/20429-38 为验收。经 15 轮（R18-32）系统诊断（含 2 子代理 + kernel/host 双探针 + 大 cap 对照），确认该缺陷复杂度超出零散轮次可安全修复范围，正式转入**专项任务清单**（不在此占用每轮推进其他交付）。
- **2026-09-03（Round 31）**：止损（guarded per-region assemble 保留、双后端绿）；①coalesce/图窗 ②per-region 预算 ③多 group 下游取数三规格归档。
- **2026-09-03（Round 30）**：修复推进与止损。落地**guarded per-region assemble**（仅 calling_regions 间存在真实空隙时按 region 满预算分区拼接；连续/合并 span 保持单次调用）——双后端 4/4（hc-broad/kmer-list/streaming-contract/mutect2-contract）回归通过。**残余两级缺陷已精确定位**：(1) 相邻 AssemblyRegion 被 coalesce/图窗合并成单个长跨度（activity 连续跨 interval 间隙时 coalesce 的 ~500bp 合并界把 200+bp 间隙也并掉）→ 单种子 + 预算跨区耗尽，B 区真实变异仍无候选（full 窗实测 candidate_sites 320 但仅 4 calls）；(2) 即便候选齐备，多 group/分区下游基因型化仍只出首尾子集（20228/20429 等 24-31 reads 支持的真实变异仍漏）。**下一步规格**：① coalesce 合并界应受 max_assembly_region_size+2*halo 约束、避免把应分隔的 region 并成长跨度；② assemble 预算改 per-calling-region（regions>1 时无条件分区）并同步更新 broad oracle 的 candidate_softclip_suppressed 语义断言（VCF 侧与 GATK 一致时以新语义为准）；③ 修复多 group 下游基因型化只取首尾 group 的问题。修复后以 fixtures/chr20 全窗（期望含 20228/20229/20429-38 等真实 calls）+ 双后端 210/210 回归 + 防复发 CTest 验收。临时探针已撤净，源码 canonical（仅含 guarded 增量），双后端绿。
- **2026-09-03（Round 29）**：诊断闭环 + 止损。仪器化（env 门控）实测结论：**凡 reference 窗口数 >1（多 interval 或多 AssemblyRegion 合并窗），HC 下游全部候选基因型化失效输出 0 calls**（连首窗口也丢），单窗口一切正常；graph 候选派生与路径枚举均非根因；`append_graph_variant_candidates` 双配置均 appended=0，calls 全来自 pileup 候选→PairHMM→基因型化路径（run() 7855-7990 的 PL/confidence 段），嫌疑集中于按候选 PL 行派生在多组/分区下的全灭。临时探针已全部撤销，源码恢复 canonical，双后端重建中。**归档为独立修复任务（A-line 阻断缺陷）**：需按 calling_region 逐窗审计 run() 候选→PairHMM→基因型链路（修复候选：逐窗口全流程或分区模式修正），完成后以 fixtures/chr20 复现 + 双后端 210/210 回归 + 防复发 CTest 验收。
- **2026-09-03（Round 28）**：仪器化：graph 候选派生双配置均 0；窗口数>1 全灭（连首窗）；排除派生/枚举/分区/合并。
- **2026-09-03（Round 27）**：探针证明相邻区间合并为 800bp 单窗单种子（DFS 深度 256 够不到下游变异）；远距双窗各自出路径（2+8）但输出 0 → 多窗口下游全丢。
- **2026-09-03（Round 26）**：Host read/ref 拼接逐行审阅正确；kernel 多段支持健全；max_paths 未触发。
- **2026-09-03（Round 25）**：路径播种遍历所有窗口；上限非因。
- **2026-09-03（Round 24）**：bug 定位推进到 **k-mer graph 多窗口路径枚举**：A+B 两窗口合并运行时 graph_nodes=1360（相加，两窗口 reads 均入图）但 graph_haplotype_paths=2（仅窗口 A）；B 单独运行 paths=6/nodes=1050。candidate_sites 跨窗口正确（64），B 区候选因无图单体型而无 PairHMM 似然 → 无 call。kernel 层（kmer_graph.cpp build_kmer_graph_kokkos）对 reference_offsets>2 的多段输入疑似只在首段做 ALT 锚收集/路径遍历；run()（calling_pipeline.cpp 7100-7240）将窗口拼接为单 graph_input。子代理 2543a9cf 携此精确线索定位行号与机制。
- **2026-09-03（Round 23）**：双 interval telemetry 决定性证据：B 完整处理（reads 151/regions 3/haplotypes 71）但 variant_calls=2 仅 A → 丢失在最终构建段；run() 位于 calling_pipeline.cpp:6644，候选 ~7236，观察切片 ~7455+。子代理 900ba714 已中断（漫游无果）。
- **2026-09-03（Round 22）**：确认丢弃在 interval/chunk 级；merge_stream_result clear payload 后不 append（hc_call.cpp:1476）；主循环 sink ~4523。calling::run 单 contig 一次调用（hc_call.cpp:4516）。
- **2026-09-03（Round 21）**：默认路径单 contig→单次 calling::run；丢弃点在 pipeline 内部。
- **2026-09-03（Round 20）**：两个分离 -L 只出首个 interval calls；telemetry 标量全量求和而 payload 仅首 chunk。
- **2026-09-03（Round 19.5）**：bug 范围收窄至 **interval/chunk 级**（非 AssemblyRegion 级）：单次运行 `-L A -L B` 两个分离 interval 只输出首个 interval 的 calls；多 region 单 interval 同理。telemetry 标量跨 chunk 求和正常（candidate_sites=64 含全部）而 payload 向量仅首 chunk（variant_calls=2）。`merge_stream_result`（hc_call.cpp:1476）首次 merge 后 **clear 全部 per-chunk payload 向量且后续 merge 从不 append**（注释"已物化进输出 body"）；主循环 sink（hc_call.cpp ~4523）对每个 encoded batch 调 vcf_text→writer.write。子代理继续定位 decode 段为何 chunk≥2 无记录。
- **2026-09-03（Round 19.5）**：合并窗口只出首个 region calls；merge_stream_result 语义为旁证。
- **2026-09-03（Round 18.5）**：确认真实数据 bug（多 region/多 interval 丢后续 calls；read 证据判真）。
- **2026-09-03（Round 18）**：**确认真实数据功能 bug**（chr20 fixture）：native HC 在 `-L 20:10019901-10020500`（跨 ≥2 AssemblyRegion）只输出首个 region 的 calls（19967/19969）；同一 reads 隔离运行 `-L 20:10020200-10020500` 能正常输出 10020228 G>C（24/24 reads 支持的真实变异）等。**多 AssemblyRegion 时后续 region 的 calls 被静默丢弃**——read 证据均为真实变异（10020228 24×C、10020429 31×G、10020434 32×T 等），与 GATK 互有错漏（GATK 在 19967/19969/20431/20434 的 ALT 无 reads 支持）。OpenMP/Serial 一致。已委派子代理在 calling_pipeline.cpp 定位丢弃点（group_ordinal 循环 ~4188 起）。
- **2026-09-03（Round 18）**：chr20 分歧诊断。验证参考序列在该区正确（ref20mnp @10019967=C，与 hg19/UCSC 及 ACGTN 流双重一致）；**位点 read 解剖**（23 条 overlap，22×T + 1×A，4 条 dup 标记）：reads 压倒性支持 **T** → native 的 ALT=T **与 read 证据一致**；GATK 的 ALT=G（QUAL 70.64、DP 4）**无任何 read 支持**——该分歧更可能是 GATK 在 10019967/69 密集 SNP 区的装配/低深度伪 call（或 DP 语义差异：GATK INFO DP=4 vs native 18，native=23-5dup 过滤口径）。**结论**：此样本不构成 native 明确 bug；作为 A1 深挖靶点保留（需 GATK 侧逐 read 装配追踪验证）。fixture 与双工具输出已留存 `fixtures/chr20/`；双后端基线不变。
- **2026-09-03（Round 17）**：**A3 关键进展：第二个真实变异 fixture 落地 + 发现真实分歧**。全量（113MB）匿名 LFS 下载 b37 chr20/21 参考成功；清洗串后以 UCSC 双锚点（1,000,000 与 10,020,000）验证坐标映射 **fs20[i]=hg19 base i**（文件首部冗余一字符导致行内偏移公式差 1，清洗串+锚点法绕开）；取真实 chr20 10,019,901–10,021,411 区序列合成全长度参考（`fixtures/chr20/ref20mnp.fasta` + dict/fai，chr21 全 N）。真实 mnp.bam（405 mapped reads，含 MNP，GATK BuildBamIndex 建索引需显式 `-O`）。**native 与 GATK 4.6.2.0 HC 在该真实区出现分歧**：chr20:10019967/10019969 等处 GATK ALT=G（DP 4）vs native ALT=T/A（DP 18）——REF 一致但 ALT/深度不同，OpenMP=Serial 复现（非随机）；native 只出 2 calls vs GATK 13（相邻密集 SNP 区）。这是 A3/A1 首个 chr17 之外的真实分歧样本，是**下一步诊断的靶点**（候选过滤/DP 差异 vs 装配差异待分）。fixture 已自洽收编：`fixtures/chr20/{ref20mnp.fasta,dict,fai, mnp.bam,bai}` + 双方 VCF。双后端基线仍 210/210。
- **2026-09-03（Round 15/16）**：**A3 真实 fixture 扩展的数据获取实证**。盘点：`large/` 目录全部为 Git-LFS 指针（CEUTrio 等 131–133B）；树内真实数据仅 chr17 69-70k 一套。验证了 **匿名 LFS media 下载可行**（`media.githubusercontent.com` HTTP 200/range 可断点），并用 HTTP Range 抓取 b37 chr20 1Mb 区真实碱基（与 UCSC hg19 API 双交叉验证坐标，锚定位点 `1,000,000=T`），合成 chr20/chr21 全长度参考（113MB、dict/fai、占位 M5 已去避免校验冲突）。**实证结论**：树内 CEUTrio b37 ch20 BAM（1m-1m1k/4379150/10m-10m100）均为 **BQSR 洁净区**（ch20 1m-1m1k 运行：858 读、graph 0 节点、0 变异）——只适合 BQSR/比对类 oracle，**不适合作为含变异调用的 HC 语料**；A3 的 ≥5 真实变异 fixture 需另取真实 WGS 变异区数据（如 NA12878 chr20 真实 reads + 位点，LFS 大文件或外部下载，网络/带宽为约束）。已记录配方（Range 抓取 + 参考合成 + 双源坐标校验），本轮无代码/测试改动，双后端基线仍 210/210。
- **2026-09-03（Round 14）**：VR 内部数值一致性收口。原 host 端 `mixture_log_probability`（bad-variant 选择与 missing-annotation 边际化用）仍是自然对数域（`ln w - 0.5(q+lndet)` 再 `/ln10`），与已转 log10 域的 score kernel 不一致；重写为 Java `evaluateDatum` 的 log10 组合（两阶段 crossProdTmp 二次型 + `log10(w)` + `(-0.5q)/log10` + `eval_log10_denominators` + `log10sumLog10`），`select_worst_variants` 的 cutoff/排序与 `marginalized_log10_probability`（`pow(10,·)` 累加）现在与 kernel 用同一数值定义。双后端 7/7 VR oracle 全过，vbem raw-bit gate（VQSLOD delta 0.0）保持。勘查结论：A3 的"≥5 真实 fixture"受参考资源限制（现有 HC fixture 仅 chr17 1Mb 一套；CEUTrio b37 chr20/21 WGS BAM 在本地但缺配套 b37 参考）；A5 contamination-joint 差 0.045 属 release-specific 校准的实质算法差异，非快赢。
- **2026-09-03（Round 13）**：A 线大任务勘查（无代码改动，双后端基线仍为 Round 12 锁定的 210/210）。flow/pairhmm 定向 CTest 复证（flow 5/5 含 4-case GATK raw-bit oracle、pairhmm 3/3）。A2 剩余主项 = GKL flow corpus raw-bit 扩展；阻塞点已定位：GATK 从零构造 flow read 的 Java harness 接近 mini flow 引擎规模（GATK 自带 `FlowBasedPairHMMUnitTest.testComputeLikelihoods` 因"difficult to generate ArtificialFlowBasedReads from scratch"而 disabled），扩展需以 sample.bam/sample.t0.bam 真实语料驱动完整 flow assembly。A1 同为大工程（GATK `ReadThreadingGraph` 内部 raw-bit oracle 尚不存在）。Track D 需 GPU/SLURM 环境（本机无 nvidia-smi/srun/nextflow）。审计文档已补 Round 5–11 汇总（见 `PROGRESS_AUDIT_2026-08-31.md` 顶部）。
- **2026-09-03（Round 11）**：**B1 完成：VQSLOD raw-bit**。E-step/score kernel 改为 Java log10 域组合：`refresh_variational_cache` 每项独立除以 log(10)（cachedDenomLog10）；responsibility kernel 与 score kernel 采用两阶段 crossProdTmp 二次型（`(x_j-μ_j)·Inv[j][i]` 转置序）+ `pow(10,·)` 归一 / `log10sumLog10`（`max + log10(1+Σ_{i≠max}10^(x_i-max))`）；新增 per-component `eval_log10_denominators` = `log10(pow(2π,-N/2)) + log10(pow(det,-0.5))`。**决定性发现**：VQSLOD 5e-5 残差的真正来源是**序列化**——GATK `VariantDataManager` 把 VQSLOD 作为 `String.format("%.4f", datum.lod)` **字符串**写入 INFO（非 float），而 native 用 float 编码 + 6 位有效数字写出；native 改为 double 保留 + `%.4f` 文本写出后 **VQSLOD delta = 0.0**。oracle 断言从 1e-3 收紧为**精确 token 相等（raw-bit through serialization）**；双后端 7/7 VR oracle（vbem/gatk-model/annotation-order/sample-every/zero-variance/attempt-iteration/contract）全过；compute_progress_score 与双后端 fixture-digest/progress-score CTest 门禁通过。模型 report 残余 5.2e-9（≤1e-6；正模型迭代路径的边际 raw-bit 化收益低，记录为已知边界）。
- **2026-09-03（Round 10）**：**B1 第二阶段：Jama LU 位级移植 + refresh_precision 接入**。确认 GATK 4.6.2.0 捆绑的是 **JAMA 1.0.3**（2012-11-14，left-looking Crout/Doolittle + LUcolj 列主元），`Matrix.inverse() = LU.solve(identity)`、`Matrix.det() = LU.det()`。交付 `fastgatk-native/src/vqsr_jama_lu.hpp`（按反编译结构移植 decompose/det/solve，含 pivot 行置换）；**随机 SPD 矩阵探针 300/300 值（逆矩阵+行列式）位级一致**。variant_recalibrator_tool.cpp：Gauss-Jordan `invert_matrix` 替换为 `invert_jama_lu`，`log_determinants` 语义从 Σlog|pivot| 改为 **log(Π pivot)（Java 先求 LU 积再 Math.log）**。结果：**负模型参数 ≤2e-15（~4 ulp）、正模型均值/协方差残余 3.6e-9/5.2e-9**（此前正模型 ~1e-9 相对误差级）；pMix 表位级一致；双后端 7/7 VR oracle 全过（OpenMP + Serial）。VQSLOD 仍 5e-5 —— 差异定位在 scoring/E-step 结构域（native 自然对数 vs Java log10/pow 组合 + 两阶段 crossProdTmp 二次型 + pMixLog10 精确归一 + log10sumLog10 + eval cachedDenom 的 2π/det 对数项）。**下一轮清单**：responsibility E-step 与 score kernel 改为 Java log10 域组合（含 eval denom = log10(pow(2π,-N/2)) + log10(pow(det,-0.5)) 与 pMixLog10 存储）。
- **2026-09-03（Round 7）**：**B1（VBEM raw-bit）第一阶段开工**。建立 raw-bit 对照方法：glibc log/log10/pow 与 HotSpot 在探针集上位级一致（NaN 符号位除外）；发现 native digamma 为近似实现（recurrence 到 8 + B8 渐近、`1.0/252.0` 乘）而 GATK 4.6.2.0 用的是 **commons-math3 3.5 的递归 digamma（C_LIMIT=49、`inv/252.0` 除、FastMath.log）**。交付物：`fastgatk-native/src/vqsr_fastmath_tables.hpp`（反射转储 FastMath.log 全部 1024×2 LN_MANT + 系数表，bit_cast 保位）+ `vqsr_fastmath.hpp`（FastMath.log 与 digamma 精确移植，**探针 17/17 位级一致**）。修复组 A：M-step sumProb 基数（Java 先 1e-10 再逐 datum 加，native team0 累加器注入 base）、协方差 regCovar 单步结合序（Java `prod+reg` 先加再累加）、final 模型 `timesEquals(1/sumProb)` 倒数乘法。结果：vbem oracle `max_abs_model_parameter_delta` ≈ **5.2e-9**（此前仅知 <1e-6），双后端 7/7 VR oracle（contract/gatk-model/vbem/annotation-order/sample-every/zero-variance/attempt-iteration）全过、无回归。**教训**：final 权重去 clamp 导致 jitter/零方差 fixture 出现 0/0→NaN→scored 0，已还原（保留 1e-10 下限，作为近空组件的数值护栏）。**遗留**：VQSLOD max delta 仍 5e-5（模型已近 raw-bit，差异在 scoring 路径：native 自然对数域 vs Java log10/pow 域评估器 + pMix 归一 + nanTolerantLog10SumLog10 + Jama LU det/inverse，下一轮按此清单推进）。
- **2026-09-03（Round 5）**：**全部 48 工具晋级完毕，registry 清零 prototype**。前置：双后端干净重建后全量 CTest 基线 OpenMP 210/210、Serial 210/210（期间修复一个真实测试 bug：`fastgatk-kernels/tests/api_smoke.cpp` 的 quality_activity 用例漏填 `indel_counts` 占位字段，导致 `{0,2,4}` 错位进 indel_counts 触发参数校验异常；另修复 `verify_mutect2_mismapping_rate_boundary.py` 的 jshell prefs 只读 HOME 问题——加 `-J-Djava.util.prefs.userRoot` 密闭目录）。随后按证据标准批量晋级 16 个核心工具（BQSR 四件套、Mutect2 链路五件、CNV 五件、VQSR 两件）：每个工具补 `compatibility_scope` + `fallback_boundaries` + 事实性晋级 note（均含双后端干净构建 CTest 证据与既有 benchmark 引用；VBEM raw-bit、MCMC raw-bit、PoN HDF5 字节级差异仍为显式 pending 边界）。`verify_dispatcher.py` 12 条断言同步（142 项检查通过）；`compute_progress_score.py` 通过，`registry_snapshot` 更新为 {48 entries, 0 prototype, 47 contract-compatible, 1 adapter, fraction 0.9792}；双后端 progress-score/fixture-digest CTest 门禁通过。全局分仍为 0.5385（专家加权，不随 registry 机械变化）。**剩余主线即 Track A/B/D 的硬 parity 任务**：A1 装配图原始位级、A2 PairHMM flow、A3/A4 后验、A5 joint 模型、A6 WGS benchmark；B1 VBEM raw-bit（~4e-5）、B2 tranche joint、B3 MCMC、B4 PoN HDF5、B5 CNV e2e；D1–D5 GPU/集群。
- **2026-09-03（Round 4/5）**：**构建目录陈旧事故 + Mutect2 六连失败的根因定位**。Round 4 评估剩余 16 个 prototype 时发现 OpenMP 批次中 6 个 Mutect2 CTest 全部 SIGSEGV。系统性排查：Serial 二进制不崩；崩溃 100% 可复现且仅在 stdout 重定向时（manifest 文件 0 字节）；gdb 初判 `ResourceSnapshot::to_json`→`json_escape` 读到 size 被指针值覆盖的字符串；ASan（-O2 与 -O3）均无任何报告；`-O2`/`-g`/`-fno-omit-frame-pointer` 任一改动都使崩溃消失。最终决定性实验：**全新干净构建（同源同 flag 同 ZEN4 arch）不崩** → 结论：主 `build/` 目录为**跨多次增量构建的陈旧对象混链**（对象时间戳横跨 08-21 至 09-02，Kokkos 核心对象为 08-22 编译，期间 arch 标志与头文件多次变更；Makefile 生成器不追踪 flag 变更，旧对象不会重编译）。HC 二进制不崩是因为其对象恰好自洽。处置：两个 canonical 构建目录（`build`、`build-serial`）整体移出为 `*-stale-suspect` 保留备查，从零全量重建两个后端，随后重跑全部 CTest 双后端作为晋级证据基线。**教训**：晋级证据只接受干净构建产物；`cmake --build` 的增量正确性不足以支撑审计级证据。同期完成：verify_all.sh 补齐 8 个此前未接入的 oracle/verify 脚本 + 4 个 benchmark（`benchmark_gather_bqsr_reports.py` 为新建，合成 GATKReport v1.1 shard p50≈28ms）。
- **2026-09-02（Round 3）**：Track C 批量晋级 13 工具（C4–C10 全部完成 + 额外 3 个合格工具）。证据：OpenMP 后端 43 项相关 CTest 全过（两批 20/20 + 23/23），Serial 后端 43/43 全过；CalculateContamination 的 verify/oracle/benchmark 此前未接入 `verify_all.sh`，本轮补齐 3 行（其 GATK oracle 5 fixture bit-exact）。registry：contract-compatible 16→**31**、prototype 31→**16**、adapter 1（fraction 0.6458）；每个晋级工具补 `compatibility_scope` + `fallback_boundaries` 与事实性 benchmark/telemetry note。`verify_dispatcher.py` 12 条状态断言同步，142 项检查通过；`compute_progress_score.py` 与双后端 progress-score/fixture-digest 门禁通过。晋级工具清单：CombineGVCFs、VariantEval、ValidateVariants、LeftAlignAndTrimVariants、ReblockGVCF、CollectAllelicCounts、DepthOfCoverage、GetPileupSummaries、GatherPileupSummaries、CalculateContamination、VariantsToTable、SelectVariants、VariantFiltration。剩余 16 个 prototype：BQSR 四件套、Mutect2 链路五件、CNV 五件、VQSR 两件、GenomicsDBImport(adapter 除外）——即核心 caller/模型工具，属于 Track A/B 的硬任务。
- **2026-09-02（Round 2）**：C1（MarkDuplicates）与 C2（SortSam）完成 registry 晋级。证据：OpenMP 后端 8/8 CTest（contract + gatk-oracle + duplicate/cli-boundary）通过、Serial 6/6 通过；`verify_mark_duplicates*.py` 三个 pinned oracle、`verify_sort_sam*.py`、`benchmark_mark_duplicates.py`（p50≈30.6ms/512 records 含 spill）、`benchmark_sort_sam.py` 全部通过；`verify_dispatcher.py` 142 项检查通过（已同步晋级断言）；`compute_progress_score.py` 通过，registry_snapshot 更新为 18 contract-compatible / 29 prototype / 1 adapter（fraction 0.375）。注意：全局分数 0.5385 是专家加权工作流评分，不随 registry 机械变化；工具晋级推动的是 long_tail_tools 工作流的实际边界。晋级时同步为两个工具补写了 compatibility_scope 与 fallback_boundaries（cloud/MD5/release-specific Picard 语义仍为显式 fallback）。

## 0. 不变约束（所有任务共同遵守）

1. **完成定义**（主计划 §11）：每个任务必须同时交付 —— 头文件/API、实现、CMake target、unit/golden test、benchmark、README 更新、manifest schema 更新、fallback 说明；只交 kernel 不算完成。
2. **停止条件**（§9）：任何工具在 输出契约 / oracle / 资源安全门 失败时，状态退回 `fallback`；禁止用"性能更快"覆盖语义失败。
3. **收尾门禁（全局串行点）**：`verify_all.sh`、`compute_progress_score.py`、`pinned_fixture_digests.sha256`、双后端 CTest（OpenMP + Serial）。并行任务完成后必须串行过这道门。
4. 禁止将不同参数的热 kernel 数字宣传为端到端加速比；benchmark 必须记录 decode/prepare/kernel/collect/write 分项与同参数 Java 对照。

## 1. Track A — HC/Mutect2 核心 parity（最高权重主线）

**前置依赖关系**：A1（Assembly graph）是 A3/A4 的根；A2、A5、A6 与 A1 并行。

| # | 任务 | 关键交付物 / 验收 | 并行性 | 状态 |
|---|---|---|---|---|
| A1 | ReadThreadingAssembler 完整 parity | k-mer graph→haplotype 遍历、pruning、cycle-retry telemetry 对齐 GATK；新增 `fastgatk-hc-graph-gatk-oracle` 双后端 raw-bit；manifest 记录 selected-k/cycle | 主线根 | ☐ |
| A2 | PairHMM flow 模式 | GKL flow corpus raw-bit；4096 matrix gate 扩展；length bucket/workspace pool 复用 | **与 A1 并行** | ◐ R13 勘查：raw-bit 内核与 4 个 pinned case 已落地（`fastgatk-kernels/tests/flow_pairhmm_gatk_oracle.cpp`，双后端过，flow 5/5 + pairhmm 3/3 CTest 绿）；length bucket/workspace pool 复用已实现；**剩余 GKL corpus 扩展的阻塞点**：GATK 侧从零构造 flow read 需近似 mini flow 引擎（GATK 自带 `FlowBasedPairHMMUnitTest.testComputeLikelihoods` 即因此 disabled），须用 sample.bam/sample.t0.bam 真实语料驱动完整 flow assembly 才能产出更多 Java 参考位 |
| A3 | HC genotyping/posterior 全对齐 | allele subsetting、AF/posterior、RCM PL 已修以外边界；HC VCF/GVCF 字段级 exact corpus 扩到 ≥5 真实 fixture | 依赖 A1 | ☐ |
| A4 | Mutect2 release-specific posterior | somatic likelihood、TLOD/AF prior calibration；`verify_mutect2_gatk_oracle.py` posterior 面 exact | 依赖 A1；**与 A5 并行** | ☐ |
| A5 | FilterMutectCalls joint filter model | 现仅 orientation/ROQ bounded；joint 输出在 ≥2 个 tumor fixture 上 exact | **与 A4 并行**（分文件） | ☐ |
| A6 | HC/Mutect2 WGS 规模 amortized benchmark | ≥NA12878 30x shard，同参数 vs Java 的 wall/RSS/IO 基线（K4/K5 层） | **独立并行** | ☐ |

## 2. Track B — 模型类工具（VQSR / CNV）

| # | 任务 | 关键交付物 / 验收 | 并行性 | 状态 |
|---|---|---|---|---|
| B1 | VariantRecalibrator VBEM raw-bit parity | 当前 VQSLOD 最大差 ~4e-5 → raw-bit；收敛、校准、provenance | **【P】** | ☑ 2026-09-03 完成（范围：pinned 2-Gaussian fixture）：R7 digamma/FastMath.log 位级移植 + M-step 结合序；R10 Jama 1.0.3 LU 位级移植（探针 300/300）；R11 E-step/score kernel log10 域组合 + VQSLOD 按 GATK `String.format("%.4f")` 文本写出 → **VQSLOD delta 0.0**，oracle 从 1e-3 收紧为精确 token 相等，双后端 7/7 VR oracle 通过。模型 report 残余 5.2e-9（≤1e-6 gate；正模型迭代路径 raw-bit 化收益低，作为已知边界记录） |
| B2 | ApplyVQSR/GatherTranches 两次 tranche 联合 | 跨 scatter 的 allele 联合状态与 tranche walk | **【P】** | ☐ |
| B3 | ModelSegments 完整 Gibbs/MCMC/RNG parity | 现仅 conditional-point；Java slice-sampler raw-bit、完整 latent-indicator MCMC、credible-interval smoothing | **【P】** | ☐ |
| B4 | CreateReadCountPanelOfNormals HDF5 byte identity + Spark/cloud metadata | v7 schema byte identity；release-specific Spark metadata | **【P】** | ☐ |
| B5 | CNV 端到端 benchmark | counts→PoN→denoise→segments→call 完整链路基线 | 依赖 B3/B4 | ☐ |

## 3. Track C — 长尾工具批次（31 prototype 中选 10）

> 每个工具统一走五门生命周期：api_cli → oracle(raw-bit) → format_sidecar → resource_io → e2e_perf。

| # | 工具 | 目标 | 并行性 | 状态 |
|---|---|---|---|---|
| C1 | MarkDuplicates | external-memory spill + metrics contract → contract-compatible | **【P】** | ☑ 2026-09-02 完成：三 pinned Picard/GATK oracle + contract 双后端通过（OpenMP 4/4、Serial 3/3），spill/resume 已落地，registry 晋级 + dispatcher 142 项检查通过 |
| C2 | SortSam | 完整 Picard sort-order / MD5 / cloud → contract-compatible | **【P】** | ☑ 2026-09-02 完成：contract + gatk-oracle + duplicate-oracle + cli-boundary 双后端通过，registry 晋级（cloud/MD5/release-specific 仍 fallback_boundaries） |
| C3 | BQSR 云端分布式边界 | resource_io 65→90（依赖 runtime staging） | **【P】** | ☐ |
| C4 | CombineGVCFs / GenomicsDB | 完整 batch/locus block；GenomicsDB native workspace 或 adapter 强化 | **【P】** | ☑ 2026-09-02 完成（CombineGVCFs 晋级；GenomicsDB 保持 adapter 决策不变） |
| C5 | VariantEval | 剩余 evaluator / graph / plugin → contract-compatible | **【P】** | ☑ 2026-09-02 完成（plugin surface 仍 fallback_boundaries） |
| C6 | SelectVariants / VariantFiltration | JEXL 全量 + concordance graph | **【P】** | ☑ 2026-09-02 完成（complex JEXL/concordance 仍 fallback_boundaries） |
| C7 | ValidateVariants | 完整 Tribble/htsjdk 字典/GVCF annotation 语义 | **【P】** | ☑ 2026-09-02 完成（完整 Tribble 语义仍 fallback_boundaries） |
| C8 | CollectAllelicCounts / DepthOfCoverage | HDF5 / cloud / report parity 边界 | **【P】** | ☑ 2026-09-02 完成（cloud/HDF5 仍 fallback_boundaries） |
| C9 | LeftAlignAndTrimVariants | symbolic / cloud / bit-identical normalization | **【P】** | ☑ 2026-09-02 完成（symbolic/bit-identical 仍 fallback_boundaries） |
| C10 | ReblockGVCF | 真实 shard 组合边界收掉（现 4/4 定向基线） | **【P】** | ☑ 2026-09-02 完成（posteriors/deletion-trimming 仍 fallback_boundaries） |

## 4. Track D — GPU 与集群（从 2% 起步）

| # | 任务 | 关键交付物 / 验收 | 并行性 | 状态 |
|---|---|---|---|---|
| D1 | PairHMM CUDA/HIP correctness backend | Strict raw-bit；NaN/underflow policy 文档 | **【P】** | ☐ |
| D2 | SW CUDA correctness + OOM → CPU fallback | 逐 CIGAR oracle；GPU OOM 可回退 | **【P】** | ☐ |
| D3 | k-mer/graph primitives GPU correctness | 按收益顺序迁移 | 依赖 D1 | ☐ |
| D4 | GPU 性能门禁 | H2D/D2H 计入后 ≥ CPU 1.5x 且输入规模过阈值才默认启用；写入门禁 schema | 依赖 D1/D2 | ☐ |
| D5 | 真实多节点 SLURM + Nextflow scatter/gather/resume | 生产证据（现仅 local/单节点 smoke） | **【P】**（与 D1–D4 并行） | ☐ |

## 5. 并行拓扑（执行视图）

```
Track A:  A1 ─→ A3 ─→ A4 ─┐
          A2 ───────────────┤（A2/A5/A6 与 A1 并行）
          A5 / A6 ──────────┘
Track B:  B1 | B2 | B3 | B4（4 路并行）→ B5（收尾）
Track C:  C1 … C10（至多 10 路并行，建议 3–4 路滚动）
Track D:  D1 ─┐
          D2 ─┼→ D3 / D4    D5（独立）
          
汇聚门（全局串行）：verify_all.sh → compute_progress_score.py → digest → CTest×2
```

四条 Track（A/B/C/D）之间零依赖、可并线推进；同 Track 内按上表依赖串并行。

## 6. 里程碑

| 里程碑 | 完成项 | 预期进度 |
|---|---|---|
| **M1** | A1、A2、B1、C1–C3、D1 | hc_mutect2_bqsr → ~65%；全局 → ~57% |
| **M2** | A3、A4、B3、C4–C7、D2、D4 | 全局 → ~61% |
| **M3** | A5、B5、C8–C10、D3、D5 | 全局 → ≥65% |

## 7. 执行须知（给执行者）

- 每个任务开工前：先读 `MODULE_IMPLEMENTATION_PLANS.md` 对应模块蓝图 + 本文件任务行。
- 永远从 pinned GATK 4.6.2.0 fixture 出发写 oracle；没 oracle 的模块默认 `prototype`。
- 长尾工具（Track C）选任务时优先选已有 verify_*.py / benchmark_*.py 脚本骨架的，能直接复用生命周期。
- Track D 的 GPU 工作只标 correctness；性能收益必须过 D4 门禁后才能宣称启用。
- 每次任务完成后更新：`progress_score.json`、`PROGRESS_AUDIT_2026-08-31.md` 对应工作流行、以及本文件状态列（☐→☑）。
