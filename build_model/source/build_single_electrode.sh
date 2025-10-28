#!/bin/bash

mkdir -p ./model
bash_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cpp_src=${bash_dir}/generate_slit.cpp
execu_bin=${bash_dir}/../../bin/generate_slit

if [[ ! -x "${execu_bin}" ]]; then
    echo "Executable file not found or not executable. Attempting to compile..."  | tee -a ./result/b_model.log >&2
    if [[ ! -f "${cpp_src}" ]]; then
        echo -e "${ERROR}Source file '${cpp_src}' not found. Cannot compile.${NC}" | tee -a ./result/b_model.log  >&2
        exit 1
    fi

    g++ "$cpp_src" -o "$execu_bin" -O3 -std=c++17
    if [[ $? -ne 0 ]]; then
        echo -e "${ERROR}Compilation of Executable file failed.${NC}" | tee -a ./result/b_model.log  >&2
        exit 1
    fi

    chmod +x "${execu_bin}"
    echo "Executable file compiled successfully." | tee -a ./result/b_model.log  >&2
fi

${execu_bin} -o ./model/single_electrode.gro
gmx editconf -f ./model/single_electrode.gro -o ./model/single_electrode.pdb