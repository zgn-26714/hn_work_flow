#!/bin/bash


rm ./result/run_md.log
set -euo pipefail

bash_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
mkdir -p result

if (( ${ismakecase} == 1 )); then
    if bash "${bash_dir}"/rerun/make_rerun_cases.sh; then
        echo -e "${GREEN}"
        echo "  ✅success! make_rerun_cases.sh completed" | tee -a ./result/run_md.log >&2
        echo "  📁The cases has been created." | tee -a ./result/run_md.log >&2
        echo -e "${NC}"
    else
        echo -e "${RED}❌ error: make_rerun_cases.sh failed${NC}" | tee -a ./result/run_md.log >&2
        exit 1
    fi
fi


export START=$rerun_start
export END=$rerun_end
if (( ${ismakejob} == 1 )); then
    if bash "${bash_dir}"/sub_job.sh; then
        echo -e "${GREEN}"
        echo "  ✅success! sub_rerun.sh completed" | tee -a ./result/run_md.log >&2
        echo "  📤The cases has been created." | tee -a ./result/run_md.log >&2
        echo -e "${NC}"
    else
        echo -e "${RED}❌ error:sub_rerun.sh failed${NC}" | tee -a ./result/run_md.log >&2
        exit 1
    fi
fi

echo -e "${GREEN}🎉🎉🎉  ALL STEPS COMPLETED SUCCESSFULLY!  🎉🎉🎉${NC}" | tee -a ./result/run_md.log >&2