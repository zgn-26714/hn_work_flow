#!/bin/bash

rm ./result/run_md.log

set -euo pipefail

bash_dir=$(builtin cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
mkdir -p result
if (( ${ismakecase} == 1 )); then
    if bash "${bash_dir}"/make_charging_cases.sh; then
        echo -e "${GREEN}"
        echo "  ✅success! make_charging_cases.sh completed" | tee -a ./result/run_md.log >&2
        echo "  📁The cases has been created." | tee -a ./result/run_md.log >&2
        echo -e "${NC}"
    else
        echo -e "${RED}❌ error: make_charging_cases.sh failed${NC}" | tee -a ./result/run_md.log >&2
        exit 1
    fi
fi

if (( ${ismakejob} == 1 )); then
    if bash "${bash_dir}"/sub_job.sh; then
        echo -e "${GREEN}"
        echo "  ✅success! sub_job.sh completed" | tee -a ./result/run_md.log >&2
        echo "  📤The cases has been created." | tee -a ./result/run_md.log >&2
        echo -e "${NC}"
    else
        echo -e "${RED}❌ error:sub_job.sh failed${NC}" | tee -a ./result/run_md.log >&2
        exit 1
    fi
fi

echo -e "${GREEN}🎉🎉🎉  ALL STEPS COMPLETED SUCCESSFULLY!  🎉🎉🎉${NC}" | tee -a ./result/run_md.log >&2