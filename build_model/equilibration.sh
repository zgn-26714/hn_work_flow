#!/bin/bash

{
echo -e "${BLUE}"
echo "**********************************************************************">&2
echo "    II.   equilibration">&2
echo "⚙️  Parameters: nsteps = ${nsteps_equ}">&2
echo "**********************************************************************">&2
echo -e "${NC}"
} | tee -a ./result/b_model.log 
if [[ -f "$GMXRC" ]]; then
    source "$GMXRC"
else
    echo -e "${ERROR}GMXRC not found: $GMXRC${NC}"  | tee -a ./result/b_model.log >&2
    exit 1
fi

echo "[step 1] clear file" | tee -a ./result/b_model.log >&2
rm -f *#

echo "[step 2] genenrate .tpr (grompp)" | tee -a ./result/b_model.log >&2
mkdir -p nvt20
cp "${MDP}".mdp ./nvt20/grompp.mdp
sed -i "s/^[[:space:]]*nsteps[[:space:]]*=.*/nsteps = ${nsteps_equ}/" ./nvt20/grompp.mdp
echo q| gmx make_ndx -f ./build/pre_eq.gro -o ./nvt20/index.ndx >> ./result/b_model.log  2>&1
gmx grompp \
    -f ./nvt20/grompp.mdp \
    -c ./build/pre_eq.gro \
    -o ./nvt20/nvt20.tpr \
    -maxwarn "${maxWarn}" \
    -p "${TOP}.top" \
    -n ./nvt20/index.ndx  >> ./result/b_model.log  2>&1 \
    || { echo -e "${ERROR} gmx grompp failed${NC}" | tee -a ./result/b_model.log >&2; exit 1; }

#  mdrun 
echo "[step 3]  mdrun" | tee -a ./result/b_model.log >&2
if [[ "$GPU" -eq 1 ]]; then
    echo -e "\tUsing GPU acceleration"  | tee -a ./result/b_model.log >&2
    mdrun_cmd=(
        gmx mdrun
        -s ./nvt20/nvt20.tpr
        -deffnm ./nvt20/nvt20
        -ntmpi 1
        -ntomp "$NPOS"
        -pme gpu
        -pmefft gpu
        -nb gpu
        -tunepme no
        -v
    )
else
    echo -e "\tUsing CPU computation"  | tee -a ./result/b_model.log >&2
    mdrun_cmd=(
        gmx mdrun
        -s ./nvt20/nvt20.tpr
        -deffnm ./nvt20/nvt20
        -ntmpi 1
        -ntomp "$NPOS"
        -v
    )
fi

echo -e "\t${mdrun_cmd[*]}" | tee -a ./result/b_model.log >&2
"${mdrun_cmd[@]}" >> ./result/b_model.log  2>&1
echo -e "${GREEN}Successfully performed a balanced simulation..${NC}" | tee -a ./result/b_model.log >&2