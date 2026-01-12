#!/bin/bash

set -euo pipefail

workdir="$(pwd)"
mkdir -p "${workdir}/slowdir"
export SLOWDIR="${workdir}/slowdir/"


SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
case "$server_machine" in
    "eninstein"|"jiaocha")
        if [ -f "generate_scanfile" ]; then
            echo "exist generate_scanfile."
        else
            g++ -o generate_scanfile ${SCRIPT_DIR}/generate_scan.cpp -O3  -std=c++17
        fi
    ;;
    "4090"|"wuchao")
        if [ -f "generate_scanfile" ]; then
            echo "exist generate_scanfile."
        else
            g++ -o generate_scanfile ${SCRIPT_DIR}/remake_scanrate.cpp -O3  -std=c++17
        fi
    ;;
esac

./generate_scanfile >> ./result/run_md.log 2>&1
