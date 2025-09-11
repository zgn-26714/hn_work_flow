#!/bin/bash
echo ""
echo "**********************************************************************">&2
echo "    II.   equilibration">&2
echo "⚙️  Parameters: nsteps = ${nsteps_equ}">&2
echo "**********************************************************************">&2

if [[ -f "$GMXRC" ]]; then
    source "$GMXRC"
else
    echo "【ERROR】GMXRC not found: $GMXRC" >&2
    exit 1
fi

echo "【step 1】 clear file">&2
rm -f *.tpr
rm -f *#
rm -rf nvt20
mkdir -p nvt20
sed -i "s/^[[:space:]]*nsteps[[:space:]]*=.*/nsteps = ${nsteps_equ}/" "$MDP".mdp

echo "【step 2】 genenrate .tpr (grompp)">&2
gmx grompp \
    -f "$MDP".mdp \
    -c pre_eq.gro \
    -o ./nvt20/nvt20.tpr \
    -maxwarn "${maxWarn}" \
    -p "${TOP}.top" \
    -n "${INDEX}.ndx" > test.log 2>&1 \
    || { echo "【ERROR】gmx grompp failed" >&2; exit 1; }

#  mdrun 
echo "【step 3】  mdrun">&2
if [[ "$GPU" -eq 1 ]]; then
    echo -e "\t【INFO】Using GPU acceleration" >&2
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
    echo -e "\t【INFO】Using CPU computation" >&2
    mdrun_cmd=(
        gmx mdrun
        -s ./nvt20/nvt20.tpr
        -deffnm ./nvt20/nvt20
        -ntmpi 1
        -ntomp "$NPOS"
        -v
    )
fi

echo -e "\t【CMD】${mdrun_cmd[*]}">&2
"${mdrun_cmd[@]}" > test.log 2>&1 