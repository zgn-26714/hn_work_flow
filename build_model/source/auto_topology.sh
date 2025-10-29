#!/bin/bash

set -euo pipefail

gro=$1
echo "[ molecules ]" >> ${TOP}.top
cat tmp >> ${TOP}.top

rm -f tmp


awk '
/^structure / { 
    struct_lines[++count] = NR
    molecules[NR] = $2
}
/^number / { number_lines[NR] = $2 }
END {
    for (i = 1; i <= count; i++) {
        line = struct_lines[i]
        molecule = molecules[line]
        for (j = line + 1; j <= NR + 10; j++) {
            if (j in number_lines) {
                print molecule, number_lines[j]
                break
            }
        }
    }
}' ${packmol}.inp| tail -n +2 | sed 's/\.pdb//' >> ${TOP}.top



# echo -e "${OK}Successfully appended molecule counts from Packmol input to topology.${NC}"