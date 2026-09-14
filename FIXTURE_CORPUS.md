# 测试夹具语料说明（tracked 入口）

本文是仓库内**受版本控制**的夹具语料总入口，回答「新机器上要跑通 oracle 需要哪些数据、从哪来」。

## 分层

| 层 | 位置 | 是否入库 | 用途 | 入口/说明 |
|---|---|---|---|---|
| 小型 pinned 夹具 | `fastgatk-native/tests/`、`fastgatk-native/tests/fixtures/` | **是** | 毫秒级契约 / oracle | `fastgatk-native/tests/pinned_fixture_digests.sha256`（14 项，由 `verify_fixture_digests.py` 门禁校验，CTest 名 `fastgatk-fixture-digest-contract`） |
| 合成 MNP smoke 语料 | `fixtures/chr20/` | 否（`.gitignore`） | HC MNP/gVCF/window-invariance 对拍 | `fixtures/MANIFEST.md`（gitignored，仅本机） |
| 真实公共基准数据 | `testdata/` | 否（`.gitignore`） | 数值/精确度验证（GIAB truth、hs37d5、BQSR known-sites、真实 reads） | `testdata/README.md` + `testdata/MANIFEST.md`（gitignored，仅本机）；抓取脚本 `fastgatk-native/scripts/fetch_real_testdata.sh`（tracked） |
| GATK 上游测试资源 | `gatk-source/src/test/resources/` | 否（`.gitignore`） | 若干 oracle 的 BAM/CRAM/FASTA 输入 | 已被 `pinned_fixture_digests.sha256` 收录（`NA12878.chr17_69k_70k.dictFix.bam` 等） |

## 新机器恢复顺序

1. `third_party/`（vendored JDK17 + pinned GATK 4.6.2.0 jar + toolchains）—— 见 `third_party/` 获取方式（构建前置）。
2. 小型 pinned 夹具：随仓库检出（`fastgatk-native/tests/`），跑 `fastgatk-fixture-digest-contract` 校验 md5。
3. `gatk-source/src/test/resources/` 中已 pin 的文件：按 `pinned_fixture_digests.sha256` 补齐。
4. `fixtures/chr20/` 与 `testdata/`：按各自 `MANIFEST.md`/`README.md` 或 `fetch_real_testdata.sh` 恢复。
5. 验证：`fastgatk-native/scripts/run_regression.sh` 默认 `FASTGATK_REQUIRE_GATK_ORACLE=1`，缺失 oracle/夹具时**响亮失败**而非静默通过。

> 注：`fixtures/`、`testdata/`、`gatk-source/`、`third_party/` 下的 `MANIFEST.md`/`README.md`
> 与被忽略的数据同处一目录，因此也**不在版本控制内**；它们是给本机浏览用的。
> 若要把大数据纳入可移植分发，需引入 git-LFS 或外部下载清单，属独立决策（未做）。
