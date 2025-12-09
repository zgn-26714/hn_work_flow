#!/bin/bash

set -euo pipefail

echo "[ molecules ]" >> ${TOP}.top


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




##### 生成bulk_top.top
echo "[ molecules ]" >> ${bulk_top}.top


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
}' ${bulk_pa}.inp| tail -n +2 | sed 's/\.pdb//' >> ${TOP}.top
# echo -e "${OK}Successfully appended molecule counts from Packmol input to topology.${NC}"