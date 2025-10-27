#!/usr/bin/bash

set -euo pipefail

bash_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

execu_bin="${bash_dir}/../bin/mod_ndx"

if [[ ! -x "${execu_bin}" ]]; then
    src_file="${bash_dir}/source/module_C++/mod_ndx.cpp"
    echo "Executable file not found or not executable. Attempting to compile..."  | tee -a ./result/analyze.log >&2
    if [[ ! -f "${src_file}" ]]; then
        echo -e "${ERROR}Source file '${src_file}' not found. Cannot compile.${NC}" | tee -a ./result/analyze.log  >&2
        exit 1
    fi

    g++ "$src_file" -o "$execu_bin" -O3 -std=c++17
    if [[ $? -ne 0 ]]; then
        echo -e "${ERROR}Compilation of Executable file failed.${NC}" | tee -a ./result/analyze.log  >&2
        exit 1
    fi

    # 确保编译后文件可执行
    chmod +x "${execu_bin}"
    echo "Executable file compiled successfully." | tee -a ./result/analyze.log  >&2
fi


cas=$1
echo "Fixing duplicate molname in case${cas}..."
cp ./case${cas}/index.ndx ./case${cas}/index.ndx.bak
${execu_bin} ./case${cas}/index.ndx ./case${cas}/index.ndx2 >> ./result/analyze.log 2>&1
rm ./case${cas}/index.ndx
mv ./case${cas}/index.ndx2 ./case${cas}/index.ndx

echo -e "${YELLOW}WARNING! Duplicate molname fixed in index.Replace the original index file with *.ndx.bak.${NC}"