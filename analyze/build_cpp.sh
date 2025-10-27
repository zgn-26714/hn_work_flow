#!/usr/bin/bash

set -euo pipefail

execu_bin=$2
src_file=$1
if [[ ! -x "${execu_bin}" ]]; then
    echo "Executable file not found or not executable. Attempting to compile..."  | tee -a ./result/b_model.log >&2
    if [[ ! -f "${src_file}" ]]; then
        echo -e "${ERROR}Source file '${src_file}' not found. Cannot compile.${NC}" | tee -a ./result/b_model.log  >&2
        exit 1
    fi

    # 尝试编译
    gmx_path=$(dirname ${build_gmx})
    installPrefix=${gmx_path}/..
    libdir=lib64
    # remove `-lonfly` if it is not onfly version
    g++ "$src_file" -o "$execu_bin" -I ${installPrefix}/include -L ${installPrefix}/${libdir} -lgromacs -O3 -std=c++17 -lonfly
    if [[ $? -ne 0 ]]; then
        echo -e "${ERROR}Compilation of Executable file failed.${NC}" | tee -a ./result/b_model.log  >&2
        exit 1
    fi

    # 确保编译后文件可执行
    chmod +x "${execu_bin}"
    echo "Executable file compiled successfully." | tee -a ./result/b_model.log  >&2
fi










