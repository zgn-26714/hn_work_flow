#!/bin/bash

set -euo pipefail

gro=$1
if echo q | gmx make_ndx -f ${gro} 2>/dev/null | grep -E "^[[:space:]]*[0-9]+" | awk '{print $2 " : " $(NF-1) " atoms"}' > index_groups.txt; then
    echo -e "${OK}Index groups generated successfully.${NC}"
else
    echo -e "${ERROR}Failed to generate index groups.${NC}"
    exit 1
fi

echo "[ molecules ]" >> ${TOP}
if grep -v -E "(System|Other)" index_groups.txt | awk -F' : ' '{print $1 " " $2}' | awk '{print $1, $2}' >> ${TOP}; then
    echo -e "${OK}Successfully updated topology with molecule counts.${NC}"
else
    echo -e "${ERROR}Failed to update topology.${NC}"
    exit 1
fi

rm -f index_groups.txt


awk '
/^structure / { struct_lines[NR] = $2 }
/^number / { number_lines[NR] = $2 }
END {
    for (line in struct_lines) {
        molecule = struct_lines[line]
        # 在当前structure行之后寻找第一个number行
        for (i = line + 1; i <= NR + 10; i++) {
            if (i in number_lines) {
                print molecule, number_lines[i]
                break
            }
        }
    }
}' ${packmol} | tail -n +2 | sed 's/\.pdb//' >> ${TOP}



# echo -e "${OK}Successfully appended molecule counts from Packmol input to topology.${NC}"