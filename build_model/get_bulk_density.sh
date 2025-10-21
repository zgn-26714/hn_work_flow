#!/bin/bash

{
echo -e "${BLUE}"
echo "**********************************************************************"
echo "    0.   Before the formal simulation"
echo "⚙️  Parameters: nsteps = ${nsteps_bulk}"
echo "**********************************************************************"
echo -e "${NC}"
} | tee -a ./result/b_model.log >&2
if [[ -f "$GMXRC" ]]; then
    source "$GMXRC"
else
    echo -e "${ERROR}GMXRC not found: $GMXRC${NC}"  | tee -a ./result/b_model.log >&2
    exit 1
fi


bash_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

bulk_name_value=$(sh "$bash_dir"/source/bulk_packmol.sh | sed -n 's/^Final structure: \(.*\)\.gro$/\1/p' )
if [ -z "$bulk_name_value" ]; then
    echo -e "${ERROR}Error: Failed to run begin_packmol.sh or no valid output${NC}" | tee -a ./result/b_model.log >&2
    exit 1
fi


echo "[step 1] clear file" | tee -a ./result/b_model.log >&2
rm -f *#

mkdir -p bulk
cp "${MDP}.mdp" ./bulk/bulk.mdp
cp "${mini_MDP}.mdp" ./bulk/mini.mdp
sed -i "s/Pcoupl[[:space:]]*=[[:space:]]*no/Pcoupl = Berendsen/" ./bulk/bulk.mdp
sed -i '/^[[:space:]]*freezegrps/s/^/;/; /^[[:space:]]*freezedim/s/^/;/' ./bulk/bulk.mdp
sed -i '/^[[:space:]]*ewald-geometry/s/^/;/' ./bulk/bulk.mdp
sed -i "s/^[[:space:]]*nsteps[[:space:]]*=.*/nsteps = ${nsteps_bulk}/" ./bulk/bulk.mdp

sed -i "s/Pcoupl[[:space:]]*=[[:space:]]*no/Pcoupl = Berendsen/" ./bulk/mini.mdp
sed -i '/^[[:space:]]*ewald-geometry/s/^/;/' ./bulk/mini.mdp
sed -i '/^[[:space:]]*freezegrps/s/^/;/; /^[[:space:]]*freezedim/s/^/;/' ./bulk/mini.mdp

if gmx grompp \
    -f "./bulk/mini.mdp" \
    -c "${bulk_name_value}.gro" \
    -p "${bulk_top}.top" \
    -o ./bulk/mini_bulk.tpr \
    -maxwarn "${maxWarn}" >> ./result/b_model.log 2>&1 ; then
    echo -e "${OK}grompp completed successfully" | tee -a ./result/b_model.log >&2
else
    echo -e "${ERROR}grompp failed.${NC}" | tee -a ./result/b_model.log >&2
    exit 1
fi


echo "[STEP 2]Running mdrun for pre-equilibration" | tee -a ./result/b_model.log >&2
if gmx mdrun \
    -s ./bulk/mini_bulk.tpr \
    -deffnm ./bulk/mini_bulk \
    -ntmpi 1 \
    -ntomp "$NPOS" \
    -v >> ./result/b_model.log 2>&1 ; then
    echo -e "${OK}min energy completed successfully" | tee -a ./result/b_model.log >&2
else
    echo -e "${ERROR}min energy failed.${NC}" | tee -a ./result/b_model.log >&2
    exit 1
fi


echo "[step 3] genenrate .tpr (grompp)" | tee -a ./result/b_model.log >&2
echo q | gmx make_ndx -f ./bulk/mini_bulk.gro -o ./bulk/bulk.index >> ./result/b_model.log  2>&1

gmx grompp \
    -f ./bulk/bulk.mdp \
    -c ./bulk/mini_bulk.gro \
    -o ./bulk/bulk.tpr \
    -maxwarn "${maxWarn}" \
    -p "${bulk_top}.top" \
    -n ./bulk/bulk.index  >> ./result/b_model.log  2>&1 \
    || { echo -e "${ERROR}gmx grompp failed${NC}" | tee -a ./result/b_model.log  >&2; exit 1; }

#  mdrun 
echo "[step 3]  mdrun" | tee -a ./result/b_model.log >&2
if [[ "$GPU" -eq 1 ]]; then
    echo -e "\tUsing GPU acceleration"  | tee -a ./result/b_model.log >&2
    mdrun_cmd=(
        gmx mdrun
        -s ./bulk/bulk.tpr
        -deffnm ./bulk/bulk
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
        -s ./bulk/bulk.tpr
        -deffnm ./bulk/bulk
        -ntmpi 1
        -ntomp "$NPOS"
        -v
    )
fi

echo -e "\t[CMD]${mdrun_cmd[*]}" | tee -a ./result/b_model.log >&2
"${mdrun_cmd[@]}" >> ./result/b_model.log  2>&1



echo "[STEP 4]Computing density profile" | tee -a ./result/b_model.log >&2
if echo 0|gmx density \
    -f ./bulk/bulk.xtc \
    -s ./bulk/bulk.tpr \
    -n ./bulk/bulk.index \
    -sl 100 \
    -b "$BE_TIME" \
    -o ./bulk/density.xvg >> ./result/b_model.log 2>&1; then
    echo -e "${OK}Density calculation completed" | tee -a ./result/b_model.log >&2
else
    echo -e "${ERROR}gmx density execution failed${NC}" | tee -a ./result/b_model.log >&2
    exit 1
fi

# =============================
# Step 4: Extract average density in Z range
# =============================
echo "[STEP 6] Analyzing density in Z range [1 - 4] nm" | tee -a ./result/b_model.log >&2
density=$(gmx analyze -f ./bulk/density.xvg -b 1 -e 4 2>/dev/null | grep "^SS1" | awk '{print $2}')

echo "OUTPUT: $density" 