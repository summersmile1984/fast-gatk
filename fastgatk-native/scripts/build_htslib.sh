#!/usr/bin/env bash
set -euo pipefail

project_root=$(cd "$(dirname "$0")/../.." && pwd)
source_root="$project_root/third_party/sources"
build_root="$project_root/third_party/htslib-build"
zlib_source="$build_root/zlib-src"
htslib_source="$build_root/htslib-src"
install_root="$build_root/install"

zlib_archive="$source_root/zlib-1.3.1.tar.gz"
htslib_archive="$source_root/htslib-1.22.1.tar.bz2"
test "$(sha256sum "$zlib_archive" | cut -d' ' -f1)" = \
     9a93b2b7dfdac77ceba5a558a580e74667dd6fede4585b91eefb60f03b72df23
test "$(sha256sum "$htslib_archive" | cut -d' ' -f1)" = \
     3dfa6eeb71db719907fe3ef7c72cb2ec9965b20b58036547c858c89b58c342f7

mkdir -p "$source_root" "$build_root" "$zlib_source" "$htslib_source"
if [ ! -f "$zlib_source/configure" ]; then
    tar -xzf "$zlib_archive" -C "$zlib_source" --strip-components=1
fi
if [ ! -f "$htslib_source/configure" ]; then
    tar -xjf "$htslib_archive" -C "$htslib_source" --strip-components=1
fi

if [ ! -f "$install_root/include/zlib.h" ]; then
    cd "$zlib_source"
    ./configure --static --prefix="$install_root"
    make -j"$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 2)"
    make install
fi

if [ ! -f "$htslib_source/libhts.a" ]; then
    cd "$htslib_source"
    CFLAGS="-I$install_root/include" LDFLAGS="-L$install_root/lib" \
      ./configure --prefix="$install_root" --disable-bz2 --disable-lzma
    make -j"$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 2)" lib-static htslib_static.mk
fi

echo "HTSlib ready: $htslib_source/libhts.a"
