#!/bin/bash

bash_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
step_cpt="${bash_dir}/.step.cpt"
result_log="./result/b_model.log"
workflow_succeeded=0

append_step_marker() {
    local marker="$1"
    echo "$marker" >> "$step_cpt"
}

detect_resume_step() {
    local step=1
    if [[ -f "$step_cpt" ]]; then
        if grep -q "^STEP_GET_FRAMES_DONE$" "$step_cpt"; then
            step=4
        elif grep -q "^STEP_EQUILIBRATION_DONE$" "$step_cpt"; then
            step=3
        elif grep -q "^STEP_BUILD_MODEL_DONE$" "$step_cpt"; then
            step=2
        fi
    fi
    echo "$step"
}

restore_original_input_files() {
    if [[ -n "${TOP:-}" && -f .top ]]; then
        cp .top "${TOP}.top"
    fi
    if [[ -n "${packmol:-}" && -f .inp ]]; then
        cp .inp "${packmol}.inp"
    fi
}

cleanup_on_exit() {
    local exit_code=$?
    if [[ $exit_code -ne 0 && $workflow_succeeded -ne 1 ]]; then
        restore_original_input_files
    fi
}

trap cleanup_on_exit EXIT

resume_step=$(detect_resume_step)
if [[ "$resume_step" -eq 1 ]]; then
    rm -rf bulk
    rm -rf build
    rm -rf nvt20
    rm -rf initial
    rm -rf result
    rm -f .top
    rm -f .inp
    mkdir -p result
    : > "$step_cpt"
else
    mkdir -p result
    echo -e "${YELLOW}[INFO]Checkpoint detected in ${step_cpt}, resume from step ${resume_step}.${NC}" | tee -a "$result_log" >&2
fi

if [[ ! -f "$GMXRC" ]]; then
    echo -e "${RED}[ERROR]GMXRC file not found: $GMXRC${NC}" | tee -a "$result_log" >&2
    exit 1
fi
source "$GMXRC"

if [[ "$resume_step" -ge 4 ]]; then
    echo -e "${GREEN}[OK]All build_model steps were already completed according to ${step_cpt}.${NC}" | tee -a "$result_log" >&2
    exit 0
fi

if [ -z "${set_density:-}" ] && [[ "$resume_step" -le 1 ]]; then
    output=$(sh "${bash_dir}/get_bulk_density.sh")

    if [ $? -ne 0 ] || ! echo "$output" | grep -q "^OUTPUT:"; then
        echo -e "${ERROR} get_bulk_density.sh failed${NC}" | tee -a "$result_log" >&2
        exit 1
    else
        echo -e "${OK}success get bulk density ${output}" | tee -a "$result_log" >&2
        density=$(echo "$output" | grep "^OUTPUT:" | cut -d' ' -f2)
        if ! [[ "$density" =~ ^[+-]?[0-9]*\.?[0-9]+([eE][+-]?[0-9]+)?$ ]]; then
            echo -e "${ERROR}Error: Invalid density value returned: $density${NC}" | tee -a "$result_log" >&2
            echo -e "${ERROR} get_bulk_density.sh failed${NC}" | tee -a "$result_log" >&2
            exit 1
        else
            echo -e "${OK}success WRITE bulk density ${density}" | tee -a "$result_log" >&2
        fi
    fi

    density=$(printf "%.10f" "$density")
    export set_density="${density}"
elif [ -z "${set_density:-}" ] && [[ "$resume_step" -gt 1 ]]; then
    echo -e "${YELLOW}[WARNING]set_density is empty, but workflow is resuming from step ${resume_step}. Continue with existing intermediate files.${NC}" | tee -a "$result_log" >&2
fi

if [[ "$resume_step" -le 1 ]]; then
    if sh "${bash_dir}/build_model.sh"; then
        append_step_marker "STEP_BUILD_MODEL_DONE"
        echo -e "${GREEN}"
        echo "##############################################" >&2
        echo "  I. success! build_model.sh completed" >&2
        echo "  The bulk density has been equilibrated." >&2
        echo "##############################################" >&2
        echo -e "${NC}"
        rm -f step*.pdb
    else
        echo -e "${ERROR} build_model.sh failed${NC}" | tee -a "$result_log" >&2
        exit 1
    fi
else
    echo -e "${YELLOW}[INFO]Skip Step I (build_model.sh), already marked done.${NC}" | tee -a "$result_log" >&2
fi

if [[ "$resume_step" -le 2 ]]; then
    if sh "${bash_dir}/equilibration.sh"; then
        append_step_marker "STEP_EQUILIBRATION_DONE"
        echo -e "${GREEN}"
        echo "##############################################" >&2
        echo "  II. success! equilibration.sh completed" >&2
        echo "##############################################" >&2
        echo -e "${NC}"
    else
        echo -e "${ERROR} equilibration.sh failed${NC}" | tee -a "$result_log" >&2
        exit 1
    fi
else
    echo -e "${YELLOW}[INFO]Skip Step II (equilibration.sh), already marked done.${NC}" | tee -a "$result_log" >&2
fi

if [[ "$resume_step" -le 3 ]]; then
    if sh "${bash_dir}/get_frames.sh"; then
        append_step_marker "STEP_GET_FRAMES_DONE"
        if [[ -f .top ]]; then
            cp .top "${TOP}.top"
        fi
        if [[ -f .inp ]]; then
            cp .inp "${packmol}.inp"
        fi
        workflow_succeeded=1
        echo ""
        echo -e "\033[1;32mALL STEPS COMPLETED SUCCESSFULLY\033[0m" >&2
        echo ""
    else
        echo -e "${ERROR} get_frames.sh failed${NC}" | tee -a "$result_log" >&2
        exit 1
    fi
else
    echo -e "${YELLOW}[INFO]Skip Step III (get_frames.sh), already marked done.${NC}" | tee -a "$result_log" >&2
fi
