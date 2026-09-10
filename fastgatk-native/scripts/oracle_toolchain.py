"""Oracle 工具链解析：统一 pinned GATK / JDK 的定位方式。

背景
----
`java` 不在 PATH 上；所有 GATK oracle 都依赖仓库内 vendored 的
`third_party/jdk17` 与 `third_party/gatk-package/gatk-4.6.2.0/`。
历史脚本各自硬编码这两个路径，其中一部分认 `JAVA` 环境变量、另一部分不认，
换机器或换 JDK 时行为不一致（见 IMPLEMENTATION_STATUS.md「工程基线」）。

约定
----
- 优先级：显式参数 > 环境变量（`JAVA` / `GATK_JAR`）> 仓库内 vendored 默认路径。
- 找不到时抛出带完整期望路径的异常，而不是让 `subprocess` 抛难懂的 OSError。
- 新写的 oracle 请使用本模块，不要再次硬编码路径。

用法
----
    from oracle_toolchain import require_toolchain
    java, gatk = require_toolchain()
    subprocess.run([java, "-Xmx1g", "-jar", str(gatk), "HaplotypeCaller", ...])
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Optional, Tuple

GATK_RELEASE = "4.6.2.0"
_REL_JDK = Path("third_party/jdk17/bin/java")
_REL_GATK = Path("third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar")
_MARKERS = ("third_party", "fastgatk-native")


def repo_root(start: Optional[Path] = None) -> Path:
    """从 start（默认本文件所在目录）向上找到仓库根。

    判定依据是同时存在 `third_party/` 与 `fastgatk-native/`，避免误停在上层目录。
    """
    here = (start or Path(__file__)).resolve()
    for candidate in [here, *here.parents]:
        if all((candidate / marker).exists() for marker in _MARKERS):
            return candidate
    raise FileNotFoundError(
        f"未能从 {here} 向上定位仓库根（需要同时存在 {' 与 '.join(_MARKERS)}）"
    )


def resolve_java(root: Optional[Path] = None) -> str:
    """返回 pinned JDK 的 java 可执行文件路径。"""
    root = root or repo_root()
    env = os.environ.get("JAVA")
    candidate = Path(env) if env else root / _REL_JDK
    if not candidate.is_file() or not os.access(candidate, os.X_OK):
        raise FileNotFoundError(
            f"找不到可执行的 java：{candidate}（可用环境变量 JAVA 覆盖；"
            f"vendored 默认路径为 {root / _REL_JDK}）"
        )
    return str(candidate)


def resolve_gatk_jar(root: Optional[Path] = None) -> Path:
    """返回 pinned GATK 4.6.2.0 的 jar 路径。"""
    root = root or repo_root()
    env = os.environ.get("GATK_JAR")
    candidate = Path(env) if env else root / _REL_GATK
    if not candidate.is_file():
        raise FileNotFoundError(
            f"找不到 pinned GATK jar：{candidate}（可用环境变量 GATK_JAR 覆盖）"
        )
    return candidate


def require_toolchain(root: Optional[Path] = None) -> Tuple[str, Path]:
    """一次性取回 (java, gatk_jar)，任何一项缺失即抛异常。"""
    root = root or repo_root()
    return resolve_java(root), resolve_gatk_jar(root)


if __name__ == "__main__":  # 便于人工自检：python3 oracle_toolchain.py
    _root = repo_root()
    _java, _gatk = require_toolchain(_root)
    print(f"repo_root : {_root}")
    print(f"java      : {_java}")
    print(f"gatk_jar  : {_gatk}")
