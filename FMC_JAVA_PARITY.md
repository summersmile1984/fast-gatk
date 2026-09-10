# FMC Java-parity 学习路径（FASTGATK_FMC_JAVA_PASSES）状态与决策档案

> 对应目标：1:1 兼容 GATK 的 C++ Host + Kokkos kernel 版本（FilterMutectCalls SomaticClusteringModel 学习链）。

## 1. 架构分工（Host + Kernel 单源数值）

- **fastgatk-kernels**（`filter_model_kokkos.{hpp,cpp}`，数值单源，全部 Kokkos 设备函数）
  - `fmc_seq_posterior_device` / `fmc_count_mixture_device`：共享 per-datum 数值原语（Java `probabilityOfSequencingError`、`logLikelihoodGivenSomatic`）。
  - `fmc_error_probabilities_kokkos`：整 pass 批量 sequencing / germline（POPAF 门控） / contamination（估计=0 退化处理）后验与组 max。
  - `fmc_model_split_quantities_kokkos`：quantile 峰分裂 / BIC 所需数量（somatic、mixture、bg softmax、quantile responsibilities）。
  - `fmc_model_em_iteration_kokkos`：单次 EM（职责快照 → 权重归一 +1 → 可选先验经验更新 → cluster 梯度/fuzzy 刷新）。
  - `gamma_math.hpp`：Commons-Math 镜像 digamma/logGamma/Dirichlet 归一（避免 host libm 漂移）。
- **fastgatk-native Host**（`fmc_java_model.hpp` + `filter_mutect_tool.cpp`）
  - 只做编排：3 次累计 pass、Java `record()` 剔除（art/ns > 0.9）、quantile 网格/峰扫描/BIC 控制流、阈值策略映射（OPTIMAL_F/CONSTANT/FDR）、结果转移进 native 模型供最终过滤。
  - std-lib 数学已全部移除（各 kernel 覆盖后删除）。

## 2. 验收面（CTest 固化，双后端可复现）

`fastgatk-filter-mutect-calls-java-*`（`FASTGATK_FMC_JAVA_PASSES=1` + pinned GATK 4.6.2.0 oracle）：

| 测试 | 验证内容 | 状态 |
|---|---|---|
| java-contamination_joint_oracle | 学习阈值/污染位点（threshold=0.1925917465371834，Δ≈0.00035 vs GATK 0.193） | pass |
| java-contamination_oracle | 污染后验/CONTQ 边界 | pass |
| java-germline_oracle | GERMQ 角点（POPAF 极端/缺失 fail-closed） | pass |
| java-normal_artifact_oracle | NALOD 边界 | pass |
| java-orientation_joint_oracle | orientation 并入学习 → threshold<1e-10 | pass |
| java-orientation_gatk_oracle | ROQ 写入 | pass |
| java-variant_index_alias_gatk_oracle | 索引别名 | pass |

kernel 级：`fastgatk-fmc-model-oracle` scenario1–5（EM、误差概率、contamination=0、分裂数量，Python a5_p2 参考一致 ~1e-13）。

## 3. 已知残差（1:1 尚未机器精确的项）

1. **contamination-joint 学习阈值 Δ≈0.00035（bg-β 分裂漂移）**：全轨迹在 EM kernel/误差概率上已验证机器级一致；残差源自分裂期 `EM(false)` 细序的 bg-β（R41/R42 定位：1.021 vs 1.045），oracle 带内（<0.05）且 Python 参考已到机器精度，C++ 侧未复现最后一位 → 待分裂/quantile 数值与 Java 逐项 diff。
2. **native 默认 mega-contract（verify_filter_mutect_calls.py）**：钉默认 INFO 驱动近似学习的遥测（observations/priors/clusters）与 AD-less INFO 重构语义 —— 属于默认路径契约，非 Java 1:1 面（R13 逐行界定：契约 1–826 行 Java 语义一致；差异从 L827 起）。
3. **真实 Mutect2 数据链**：matched Mutect2 stats + annotation 面（A5 真实 fixture）尚未建立。

## 4. 默认化决策矩阵（R15 分析）

默认路径目前走 native 近似学习（INFO 驱动）；Java-exact 学习需 `FASTGATK_FMC_JAVA_PASSES=1`。
flip 默认所需前置：
- [ ] native 默认 mega-contract 豁免/重构（差异从 L827 起，需新契约或子集划分）
- [ ] 默认 FMC 8 fixture 在 Java 语义下逐一重验 GATK pin（7 env java oracle 已覆盖其中六类语义）
- [ ] 真实 Mutect2 stats/annotation 链对齐
- [ ] 决定背驰处理：当输入缺失 AD（仅 INFO）时 Java 学习不可用（obs 语义差异）

**决策（R16 用户拍板）：维持 env 门控，Java 验收面即 1:1 交付。**
- 默认产品路径不动（native 近似学习，mega-contract 全绿 221/221 双后端）。
- `FASTGATK_FMC_JAVA_PASSES=1` + `fastgatk-filter-mutect-calls-java-*` CTest = 官方 1:1 兼容面。
- 上表前置条件不再作为默认化必须项；默认化仅当未来出现新的产品需求时重开。
- 残差（bg-β 分裂漂移 Δ0.00035 等）作为已记录边界；可选后续研究项。

## 5. 相关轮次索引
R24–R34（P1/P2 前史）· R40–R45（Java 模型 C++/env 管道）· R47/48（默认零回归、方案 B 下沉）· 新目标 R1–R15（EM/误差概率/分裂 kernel 化、bug 修复、CTest 固化、契约边界界定）。
