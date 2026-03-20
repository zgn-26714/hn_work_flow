#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"

bash "$SCRIPT_DIR"/bin/build_cpp.sh "$SCRIPT_DIR"/check_initial_files.cpp "$SCRIPT_DIR"/bin/check_initial_files

if "$SCRIPT_DIR"/bin/check_initial_files $1; then
    echo -e "${GREEN}All required files are present.${NC}"
else
    echo -e "${RED}Some required files are missing. Please check setting.log.${NC}"
    exit 1
fi
