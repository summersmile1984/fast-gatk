# Fast-GATK 会话交接状态（更新 R26 · Track A/B + 1:1 架构目标）

> 面向后续执行会话的权威状态快照；执行日志见 `A_B_EXECUTION.md`（Track A/B R1–R23 前史 + 1:1 架构目标 R1–R26）、`FMC_JAVA_PARITY.md`（FMC Java-parity 决策档案）、`NEXT_PHASE_TASKS.md`。验收口径：**本地可达 Oracle 框架**——可达项以 pinned GATK 4.6.2.0 oracle 双后端（OpenMP/Serial）机器级验证；不可达项记录为计划边界。

## 0. 当前基线（R26，全部双后端绿）
- **全量默认 CTest：OpenMP 225/225、Serial 225/225**（R24 认证树；R25/R26 零 prod 改动，宽面健康子集 31/31 复核通过）
- 新增（相对旧 213 基线）：fmc-model-oracle、fmc-trajectory-oracle(+OMP 确定性 ×3)、java-* 7 项 env GATK oracle、tlod/somatic kernels oracles 等

## 1. 已达成（机器级/常驻门禁）
- **FMC Java-parity 学习链（env 门控 `FASTGATK_FMC_JAVA_PASSES=1`）**：C++ Host 编排（3-pass、record 剔除、split/BIC 控制流、阈值策略）+ **fastgatk-kernels 单源数值**（error-probabilities、split-quantities、EM、共享 seq/mixture helper、gamma_math）。contamination-joint 学习阈值 **0.19259174653718336 vs GATK/EngineProbe 0.1925917465406033（Δ≈3.4e-12）**；R17 根因修复 = Java 每 split 跑 5×EM(false)（原 1×）。三层常驻门禁：
  1. `fastgatk-fmc-trajectory-oracle`（无 GATK 依赖，钉 bg β 1.0454669114043891/threshold 0.19259174653718336）
  2. `fastgatk-fmc-model-oracle` scenario1–5（kernel 级，Python a5_p2 参考一致 ~1e-13）
  3. `fastgatk-filter-mutect-calls-java-*`（7 项 pinned GATK oracle；contamination-joint 紧带 <1e-8；germline 内部值精确相等）
  - OMP 线程确定性：trajectory 在 1–16 线程逐位一致（CTest pin）
  - 决策（R16 用户拍板）：**维持 env 门控，Java 验收面即 1:1 交付**；默认路径保持 native 语义
- **A4 Mutect2 TLOD 引擎 parity（机器级）**：Python 转写 Java SomaticLikelihoodsEngine vs native 引擎，**max Δ=1.4e-14**（R24：先前 4.3e-7 为 VCF 8-sig 序列化测量伪影；新增 `FASTGATK_TLOD_FULL` 全精度通道 + verifier 收紧 guard <1e-8）；真实 Δ3.5@69368 = 装配/haplotype 输入边界（A1）
- **A2 可达核心**（flow PairHMM raw-bit 4-case + codec oracle）、**B1 VBEM raw-bit**、**B2 scatter→gather→apply 联合回归**：维持绿
- 204 个 `verify_*.py` GATK oracle 覆盖各工具（多数容差已紧，如 somatic-likelihood 2e-10）

## 2. 已知限制与边界（明示、勿过度宣称）
- **A-line / A1（装配/haplotype 集）**：A4 Δ3.5@69368 与多 region 连续窗缺陷根因（R15–R24 判别链收敛）；需 GATK 引擎插桩或 A1 装配对齐
- **A5 默认路径**（无 env）：native INFO 驱动近似学习 contamination-joint Δ0.042（允许带 <0.05）；native 默认 mega-contract（verify_filter_mutect_calls.py）为默认语义契约（R13 逐行界定：与 Java 路径差异自 L827 起）
- **注解层**（AS_SB_TABLE/MBQ/MMQ/MPOS 等）：语义微差（非快赢）
- 其余边界：A1/A3/A6、B3（ModelSegments MCMC/Gibbs）、B4（PoN HDF5）、B5（CNV e2e）需大型 harness/数据/长会话；Track D 需 GPU/SLURM

## 3. 资产清单（近期新增）
- `fastgatk-native/tests/fmc_trajectory_oracle.cpp` + CTest（含 `-omp1/4/16` 确定性 pin）
- `fastgatk-kernels/tests/fmc_model_oracle.cpp`（scenario1–5）
- CMake：`fastgatk-filter-mutect-calls-java-*`（7 项，FASTGATK_FMC_JAVA_PASSES=1 内建）
- `verify_mutect2_tlod_formula.py`（full-precision 通道、per_locus_deltas）、mutect2 `[FASTGATK_TLOD_FULL]` echo（FASTGATK_DEBUG_TLOD=1）
- 诊断：`FASTGATK_DEBUG_FMC` 门控 pass/EM 状态 dump（%.17g）、`.ab4/` scratch（a5_p1_full/a5_p2_reference/germ_probe/fix2 等）
- 文档：`FMC_JAVA_PARITY.md`、`A_B_EXECUTION.md`（R26）

## 4. 建议下一步（按价值；R25/R26 用户选择维持稳态）
1. 按需求指定新工具轴（方法学：kernels 单源数值 + host 编排 + 机器级 GATK oracle + 确定性 pin）
2. 长线候选（均已列明代价）：A4 multiallelic TLOD 公式面；A4 全精度通道扩展到 NLOD/后验/per-read；A1 装配 parity（整轮级）；A5 默认化前置（契约重构/豁免）
3. 回归纪律不变：每 prod 改动 → 双后端全量 + pinned oracle

## 5. 纪律（沿用不变约束）
- 晋级/宣称一律要双后端干净构建 + pinned GATK 4.6.2.0 oracle 证据；fallback 边界显式
- 不回退既有 225/225 语义；改动前先跑子集回归；引擎/注解层改动前先过 tlod-formula 与 mutect2/filter 测试
