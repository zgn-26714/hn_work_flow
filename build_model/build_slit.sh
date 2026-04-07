#!/bin/bash

# set -euo pipefail

bash_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
step_cpt="${bash_dir}/.build_slit.step.cpt"
result_log="./result/b_model.log"

append_step_marker() {
    local marker="$1"
    echo "$marker" >> "$step_cpt"
}

detect_resume_step() {
    local step=1

    if [[ -f "$step_cpt" ]]; then
        if grep -q "^STEP_GET_FRAMES_DONE$" "$step_cpt"; then
            step=6
        elif grep -q "^STEP_EQUILIBRATION_DONE$" "$step_cpt"; then
            step=5
        elif grep -q "^STEP_REPLICATE_TRANSLATE_DONE$" "$step_cpt"; then
            step=4
        elif grep -q "^STEP_BUILD_MODEL_DONE$" "$step_cpt"; then
            step=3
        elif grep -q "^STEP_SINGLE_ELECTRODE_DONE$" "$step_cpt"; then
            step=2
        fi
    fi

    echo "$step"
}

prepare_mdp_for_equilibration() {
    local bak_mdp="bak_${MDP}.mdp"

    if [[ ! -f "$bak_mdp" ]]; then
        cp "${MDP}.mdp" "$bak_mdp"
    fi

    sed -i -e 's/^freezegrps.*$/freezegrps = EL  ER  GRA/' -e 's/^freezedim.*$/freezedim   = Y Y Y Y Y Y Y Y Y/' "${MDP}.mdp"

    if [[ -f ./build/mini_re.gro ]]; then
        mv ./build/mini_re.gro ./build/pre_eq.gro
    fi
}

cleanup_temp_files() {
    rm -f ./step*.pdb
    rm -f ./bulk/*#
    rm -f ./initial/*#
    rm -f ./build/*#
    rm -f ./model/*#
    rm -f ./nvt20/*#
    rm -f bak_*
}

resume_step=$(detect_resume_step)
if [[ "$resume_step" -eq 1 ]]; then
    rm -rf model
    rm -rf bulk
    rm -rf build
    rm -rf nvt20
    rm -rf initial
    rm -rf result
    mkdir -p result
    : > "$step_cpt"
else
    mkdir -p result
    echo -e "${YELLOW}[INFO]Checkpoint detected in ${step_cpt}, resume from step ${resume_step}.${NC}" | tee -a "$result_log" >&2
fi

if [[ -f "$GMXRC" ]]; then
    source "$GMXRC" >> "$result_log" 2>&1
else
    echo -e "${ERROR}GMXRC not found: $GMXRC${NC}" | tee -a "$result_log" >&2
    exit 1
fi

if [[ "$resume_step" -ge 6 ]]; then
    echo -e "${GREEN}[OK]All build_slit steps were already completed according to ${step_cpt}.${NC}" | tee -a "$result_log" >&2
    exit 0
fi

if [[ "$resume_step" -le 1 ]]; then
    if bash "${bash_dir}/source/build_single_electrode.sh"; then
        if [[ ! -f ./model/single_electrode.pdb ]]; then
            echo -e "${ERROR} ./model/single_electrode.pdb not found${NC}" | tee -a "$result_log" >&2
            exit 1
        fi

        cp ./model/single_electrode.pdb .
        append_step_marker "STEP_SINGLE_ELECTRODE_DONE"
        echo -e "${GREEN}"
        echo "##################################################" >&2
        echo "  I. success! build_single_electrode.sh completed" >&2
        echo "##################################################" >&2
        echo -e "${NC}"
    else
        echo -e "${ERROR} build_single_electrode.sh failed${NC}" | tee -a "$result_log" >&2
        exit 1
    fi
else
    echo -e "${YELLOW}[INFO]Skip Step I (build_single_electrode.sh), already marked done.${NC}" | tee -a "$result_log" >&2
fi

if [ -z "${set_density:-}" ] && [[ "$resume_step" -le 2 ]]; then
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
elif [ -z "${set_density:-}" ] && [[ "$resume_step" -gt 2 ]]; then
    echo -e "${YELLOW}[WARNING]set_density is empty, but workflow is resuming from step ${resume_step}. Continue with existing intermediate files.${NC}" | tee -a "$result_log" >&2
fi

if [[ "$resume_step" -le 2 ]]; then
    if bash "${bash_dir}/build_model.sh"; then
        append_step_marker "STEP_BUILD_MODEL_DONE"
        echo -e "${GREEN}"
        echo "##############################################" >&2
        echo "  II. success! build_model.sh completed" >&2
        echo "  The bulk density has been equilibrated." >&2
        echo "##############################################" >&2
        echo -e "${NC}"
    else
        echo -e "${ERROR} build_model.sh failed${NC}" | tee -a "$result_log" >&2
        exit 1
    fi
else
    echo -e "${YELLOW}[INFO]Skip Step II (build_model.sh), already marked done.${NC}" | tee -a "$result_log" >&2
fi

if [[ "$resume_step" -le 3 ]]; then
    if bash "${bash_dir}/source/replicate_translate.sh"; then
        append_step_marker "STEP_REPLICATE_TRANSLATE_DONE"
        echo -e "${GREEN}"
        echo "##############################################" >&2
        echo "  III. success! replicate_translate.sh completed" >&2
        echo "##############################################" >&2
        echo -e "${NC}"
    else
        echo -e "${ERROR} replicate_translate.sh failed${NC}" | tee -a "$result_log" >&2
        exit 1
    fi
else
    echo -e "${YELLOW}[INFO]Skip Step III (replicate_translate.sh), already marked done.${NC}" | tee -a "$result_log" >&2
fi

if [[ "$resume_step" -le 4 ]]; then
    prepare_mdp_for_equilibration

    if bash "${bash_dir}/equilibration.sh"; then
        append_step_marker "STEP_EQUILIBRATION_DONE"
        echo -e "${GREEN}"
        echo "##############################################" >&2
        echo "  IV. success!  equilibration.sh completed" >&2
        echo "##############################################" >&2
        echo -e "${NC}"
    else
        echo -e "${ERROR} equilibration.sh failed${NC}" | tee -a "$result_log" >&2
        exit 1
    fi
else
    echo -e "${YELLOW}[INFO]Skip Step IV (equilibration.sh), already marked done.${NC}" | tee -a "$result_log" >&2
fi

if [[ "$resume_step" -le 5 ]]; then
    if bash "${bash_dir}/get_frames.sh"; then
        append_step_marker "STEP_GET_FRAMES_DONE"
        cleanup_temp_files
        echo -e "${GREEN}"
        echo "##############################################" >&2
        echo "  V. success!  get_frames.sh completed" >&2
        echo "##############################################" >&2
        echo -e "${NC}"
        echo ""
        echo -e "\033[1;32mALL STEPS COMPLETED SUCCESSFULLY\033[0m" >&2
        echo ""
    else
        echo -e "${ERROR} get_frames.sh failed${NC}" | tee -a "$result_log" >&2
        exit 1
    fi
else
    echo -e "${YELLOW}[INFO]Skip Step V (get_frames.sh), already marked done.${NC}" | tee -a "$result_log" >&2
fi
