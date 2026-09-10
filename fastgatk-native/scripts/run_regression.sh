#!/usr/bin/env bash
# ============================================================================
# fast-gatk 双后端回归入口
#
# 目的：把「双后端 + pinned GATK oracle + 可粘贴证据块」固化成一条命令，
#       取代手工 ctest 与 IMPLEMENTATION_STATUS.md 里不可验证的
#       “本轮改动前通过”式基线。
#
# 用法：
#   fastgatk-native/scripts/run_regression.sh                    # 双后端全量
#   fastgatk-native/scripts/run_regression.sh -R 'fastgatk-(hc|mutect2)'
#   fastgatk-native/scripts/run_regression.sh --backend omp -j 12
#   fastgatk-native/scripts/run_regression.sh --label "given-alleles 修复后"
#
# 产物：<outdir>/<backend>.log、<outdir>/evidence.md、<outdir>/summary.txt
# 退出码：0 = 全部通过；1 = 存在失败；2 = 前置条件不满足（工具链/二进制缺失）
# ============================================================================
set -uo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"

BUILD_OMP="${ROOT}/fastgatk-native/build"
BUILD_SERIAL="${ROOT}/fastgatk-native/build-serial"
CTEST="${ROOT}/third_party/toolchains/cmake-4.3.4-linux-x86_64/bin/ctest"
JAVA="${ROOT}/third_party/jdk17/bin/java"
GATK_JAR="${ROOT}/third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"

BACKEND="both"
JOBS=8
FILTER=""
LABEL=""
OUTDIR=""
TIMEOUT=""

usage() { sed -n '2,20p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    -R|--filter)  FILTER="${2:-}"; shift 2 ;;
    -j|--jobs)    JOBS="${2:-8}"; shift 2 ;;
    --backend)    BACKEND="${2:-both}"; shift 2 ;;
    --label)      LABEL="${2:-}"; shift 2 ;;
    --outdir)     OUTDIR="${2:-}"; shift 2 ;;
    --timeout)    TIMEOUT="${2:-}"; shift 2 ;;
    -h|--help)    usage ;;
    *) echo "未知参数: $1" >&2; exit 2 ;;
  esac
done

case "${BACKEND}" in omp|serial|both) ;; *) echo "--backend 必须是 omp|serial|both" >&2; exit 2 ;; esac

# ---- 工具链前置检查 --------------------------------------------------------
fail_precheck=0
[[ -x "${CTEST}" ]] || { CTEST="$(command -v ctest || true)"; }
[[ -x "${CTEST}" ]] || { echo "缺少 ctest（third_party 与 PATH 均未找到）" >&2; fail_precheck=1; }
[[ -x "${JAVA}" ]] || echo "警告：未找到 ${JAVA}（GATK oracle 可能无法运行）" >&2
[[ -f "${GATK_JAR}" ]] || echo "警告：未找到 pinned GATK jar ${GATK_JAR}" >&2

for pair in "omp:${BUILD_OMP}" "serial:${BUILD_SERIAL}"; do
  name="${pair%%:*}"; dir="${pair#*:}"
  case "${BACKEND}" in
    both) ;;
    "${name}") ;;
    *) continue ;;
  esac
  [[ -d "${dir}" ]] || { echo "缺少构建目录 ${dir}" >&2; fail_precheck=1; continue; }
  # 陈旧性检查：源码比二进制新则告警（回归结论不可信）
  newer="$(find "${ROOT}/fastgatk-native/src" "${ROOT}/fastgatk-native/include" \
                "${ROOT}/fastgatk-kernels/src" "${ROOT}/fastgatk-kernels/include" \
                -type f \( -name '*.cpp' -o -name '*.hpp' \) \
                -newer "${dir}/fastgatk-hc-call" -print -quit 2>/dev/null || true)"
  [[ -n "${newer}" ]] && echo "警告：[${name}] 源码比二进制新（${newer}），结论可能失效，请先重新构建" >&2
done
[[ "${fail_precheck}" -eq 0 ]] || exit 2

# ---- 输出目录 --------------------------------------------------------------
STAMP="$(date +%Y%m%d-%H%M%S)"
[[ -n "${OUTDIR}" ]] || OUTDIR="${ROOT}/.diag/regression/${STAMP}"
mkdir -p "${OUTDIR}"

GIT_REV="$(cd "${ROOT}" && git rev-parse --short HEAD 2>/dev/null || echo 'no-vcs')"
GIT_DIRTY="$(cd "${ROOT}" && git status --porcelain 2>/dev/null | wc -l | tr -d ' ')"

run_one() {  # $1=name $2=builddir
  local name="$1" dir="$2" log="${OUTDIR}/$1.log"
  local -a args=(--test-dir "${dir}" -j "${JOBS}" --output-on-failure)
  [[ -n "${FILTER}" ]] && args+=(-R "${FILTER}")
  [[ -n "${TIMEOUT}" ]] && args+=(--timeout "${TIMEOUT}")
  if "${CTEST}" "${args[@]}" > "${log}" 2>&1; then
    echo "PASS" > "${OUTDIR}/$1.status"
  else
    echo "FAIL" > "${OUTDIR}/$1.status"
  fi
}

# 两个后端互不共享状态（oracle 全部使用 TemporaryDirectory），可安全并行
declare -a pids=()
want() { case "${BACKEND}" in both) return 0 ;; "$1") return 0 ;; *) return 1 ;; esac; }
if want omp;    then run_one omp    "${BUILD_OMP}"    & pids+=($!); fi
if want serial; then run_one serial "${BUILD_SERIAL}" & pids+=($!); fi
for p in "${pids[@]:-}"; do [[ -n "${p}" ]] && wait "${p}"; done

# ---- 汇总 + 证据块 ---------------------------------------------------------
overall=0
{
  echo "# 双后端回归证据"
  echo
  echo "- 时间：$(date '+%Y-%m-%d %H:%M:%S %Z')"
  [[ -n "${LABEL}" ]] && echo "- 标签：${LABEL}"
  echo "- git：\`${GIT_REV}\`（未提交变更 ${GIT_DIRTY} 项）"
  echo "- 过滤：\`${FILTER:-<全量>}\`　并行度：${JOBS}"
  echo "- ctest：\`${CTEST}\`"
  echo "- pinned GATK 4.6.2.0：\`third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar\`　JDK：\`third_party/jdk17\`"
  echo
  echo "| 后端 | 构建目录 | 结果 | 通过/总数 | 耗时 |"
  echo "| --- | --- | --- | --- | --- |"
} > "${OUTDIR}/evidence.md"

for name in omp serial; do
  want "${name}" || continue
  log="${OUTDIR}/${name}.log"
  status="$(cat "${OUTDIR}/${name}.status" 2>/dev/null || echo UNKNOWN)"
  [[ "${status}" == "PASS" ]] || overall=1
  counts="$(grep -oE '[0-9]+% tests passed, [0-9]+ tests failed out of [0-9]+' "${log}" | tail -1)"
  failed="$(echo "${counts}" | grep -oE '[0-9]+ tests failed' | grep -oE '[0-9]+')"
  total="$(echo "${counts}" | grep -oE 'out of [0-9]+' | grep -oE '[0-9]+')"
  if [[ -n "${failed}" && -n "${total}" ]]; then passed=$(( total - failed )); else passed=""; fi
  # 过滤未命中任何测试时 ctest 仍返回 0，必须显式判为「无证据」而不是通过
  if [[ -z "${total}" || "${total}" == "0" ]]; then
    status="EMPTY"; overall=1
    echo "警告：[${name}] 过滤器未命中任何测试（\`${FILTER:-<全量>}\`），判为无证据" >&2
  fi
  elapsed="$(grep -oE 'Total Test time \(real\) = .*' "${log}" | tail -1 | sed 's/.*= //')"
  dir_label="$([[ "${name}" == omp ]] && echo 'OpenMP (fastgatk-native/build)' || echo 'Serial (fastgatk-native/build-serial)')"
  echo "| ${name} | \`${dir_label}\` | $([[ "${status}" == PASS ]] && echo '通过' || { [[ "${status}" == EMPTY ]] && echo '**无证据**' || echo '**失败**'; }) | ${passed:-?}/${total:-?} | ${elapsed:-?} |" >> "${OUTDIR}/evidence.md"
done

{
  echo
  echo "日志：\`$(basename "${OUTDIR}")/{omp,serial}.log\`"
  echo
  echo "> 说明：本块由 \`fastgatk-native/scripts/run_regression.sh\` 自动生成，可直接粘贴进"
  echo "> \`IMPLEMENTATION_STATUS.md\` 的「最近验证状态」。结论只对本块记录的 git 版本与二进制有效。"
} >> "${OUTDIR}/evidence.md"

cp "${OUTDIR}/evidence.md" "${OUTDIR}/summary.txt"
echo
cat "${OUTDIR}/evidence.md"
echo "证据目录：${OUTDIR}"
exit "${overall}"
