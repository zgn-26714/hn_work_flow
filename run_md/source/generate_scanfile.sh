#!/bin/bash

set -euo pipefail

workdir="$(pwd)"
mkdir -p "${workdir}/slowdir"
export SLOWDIR="${workdir}/slowdir/"


SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
g++ -o generate_scanfile ${SCRIPT_DIR}/generate_scan.cpp -O3 
./generate_scanfile >> ./result/run_md.log 2>&1
