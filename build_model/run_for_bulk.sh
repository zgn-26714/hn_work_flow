#!/bin/bash

bash_dir=$(builtin cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
step_cpt="./.run_for_bulk.step.cpt"
result_log="./result/b_model.log"

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
        elif grep -q "^STEP_PREPARE_BULK_DONE$" "$step_cpt"; then
            step=2
        fi
    fi

    echo "$step"
}

prepare_grompp_for_equilibration() {
    if [[ ! -f grompp_bak.mdp ]]; then
        cp grompp.mdp grompp_bak.mdp
    fi

    sed -i "s/Pcoupl[[:space:]]*=[[:space:]]*no/Pcoupl = Berendsen/" grompp.mdp
    sed -i '/^[[:space:]]*freezegrps/s/^/;/; /^[[:space:]]*freezedim/s/^/;/' grompp.mdp
    sed -i '/^[[:space:]]*ewald-geometry/s/^/;/' grompp.mdp
}

restore_grompp_after_completion() {
    if [[ -f grompp_bak.mdp ]]; then
        mv grompp_bak.mdp grompp.mdp
    fi
}

resume_step=$(detect_resume_step)
if [[ "$resume_step" -eq 1 ]]; then
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

if [[ ! -f "$GMXRC" ]]; then
    echo -e "${RED}[ERROR]GMXRC file not found: $GMXRC${NC}" | tee -a "$result_log" >&2
    exit 1
fi
source "$GMXRC"

if [[ "$resume_step" -ge 4 ]]; then
    echo -e "${GREEN}[OK]All run_for_bulk steps were already completed according to ${step_cpt}.${NC}" | tee -a "$result_log" >&2
    exit 0
fi

if [[ "$resume_step" -le 1 ]]; then
    output=$(sh "${bash_dir}/get_bulk_density.sh")

    if [ $? -ne 0 ] || ! echo "$output" | grep -q "^OUTPUT:"; then
        echo -e "${ERROR} get_bulk_density.sh failed${NC}" | tee -a "$result_log" >&2
        exit 1
    else
        echo -e "${OK}success get a bulk box for simulation." | tee -a "$result_log" >&2
    fi

    rm -f step*.pdb

    if [[ ! -f bulk/bulk.gro ]]; then
        echo -e "${ERROR} bulk/bulk.gro not found after get_bulk_density.sh${NC}" | tee -a "$result_log" >&2
        exit 1
    fi

    mkdir -p build
    mv bulk/bulk.gro build/pre_eq.gro
    append_step_marker "STEP_PREPARE_BULK_DONE"
    echo -e "${GREEN}"
    echo "##############################################" >&2
    echo "  I. success! " >&2
    echo "  The bulk has been pre-equilibrated." >&2
    echo "##############################################" >&2
    echo -e "${NC}"
else
    echo -e "${YELLOW}[INFO]Skip Step I (prepare bulk), already marked done.${NC}" | tee -a "$result_log" >&2
fi

if [[ "$resume_step" -le 2 ]]; then
    prepare_grompp_for_equilibration

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
        restore_grompp_after_completion
        append_step_marker "STEP_GET_FRAMES_DONE"
        : > "$step_cpt"
        echo ""
        echo -e "\033[1;32mALL STEPS COMPLETED SUCCESSFULLY\033[0m" >&2
        echo ""
    else
        echo -e "${ERROR} get_frame.sh failed${NC}" | tee -a "$result_log" >&2
        exit 1
    fi
else
    echo -e "${YELLOW}[INFO]Skip Step III (get_frames.sh), already marked done.${NC}" | tee -a "$result_log" >&2
fi
