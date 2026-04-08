#!/bin/bash
set -euo pipefail

validate_inputs() {
    if [[ ! -f ./model/single_electrode.gro ]]; then
        echo -e "${ERROR}Missing ./model/single_electrode.gro before generating Packmol input.${NC}" | tee -a ./result/b_model.log >&2
        exit 1
    fi

    read -r -a MOL_arr <<< "${MOL:-}"
    read -r -a MOLnum_arr <<< "${MOL_num:-}"
    read -r -a MOLnum_arr_bulk <<< "${MOL_num_bulk:-}"

    if [[ ${#MOL_arr[@]} -eq 0 ]]; then
        echo -e "${ERROR}MOL is empty. Cannot auto-generate Packmol input.${NC}" | tee -a ./result/b_model.log >&2
        exit 1
    fi

    if [[ ${#MOL_arr[@]} -ne ${#MOLnum_arr[@]} ]]; then
        echo -e "${ERROR}MOL and MOL_num have different lengths.${NC}" | tee -a ./result/b_model.log >&2
        exit 1
    fi

    if [[ ${#MOL_arr[@]} -ne ${#MOLnum_arr_bulk[@]} ]]; then
        echo -e "${ERROR}MOL and MOL_num_bulk have different lengths.${NC}" | tee -a ./result/b_model.log >&2
        exit 1
    fi

    for mol in "${MOL_arr[@]}"; do
        if [[ ! -f "${mol}.pdb" ]]; then
            echo -e "${ERROR}Missing molecule structure file: ${mol}.pdb${NC}" | tee -a ./result/b_model.log >&2
            exit 1
        fi
    done
}

write_packmol_header() {
    local file="$1"
    local tolerance="$2"
    local output_name="$3"

    cat > "$file" <<EOF
tolerance ${tolerance}

filetype pdb

output ${output_name}

EOF
}

append_slit_structure_block() {
    local file="$1"
    local molecule="$2"
    local count="$3"
    local zmin="$4"
    local zmax="$5"

    cat >> "$file" <<EOF
structure ${molecule}.pdb
number ${count}
inside box 1. 1. ${zmin}  ${BOX_X}  ${BOX_Y}  ${zmax}
end structure

EOF
}

append_bulk_structure_block() {
    local file="$1"
    local molecule="$2"
    local count="$3"

    cat >> "$file" <<EOF
structure ${molecule}.pdb
number ${count}
inside box 0. 0. 0. ${bulk_xbox}  ${bulk_ybox}  ${bulk_zbox}
end structure

EOF
}

validate_inputs

read -r -a MOL_arr <<< "${MOL}"
read -r -a MOLnum_arr <<< "${MOL_num}"
read -r -a MOLnum_arr_bulk <<< "${MOL_num_bulk}"

last_line=$(tail -n 1 ./model/single_electrode.gro)
read -r BOX_X BOX_Y BOX_Z <<< "$last_line"

echo "MOL is ${MOL_arr[*]}" | tee -a ./result/b_model.log >&2
echo "MOLnum is ${MOLnum_arr[*]}" | tee -a ./result/b_model.log >&2
echo "Box dimensions are: $BOX_X, $BOX_Y, $BOX_Z" | tee -a ./result/b_model.log >&2

max_ele=$(printf "%.3f" "$(bc -l <<< "$max_ele * 10")")
BOX_X=$(printf "%.3f" "$(bc -l <<< "$BOX_X * 10")")
BOX_Y=$(printf "%.3f" "$(bc -l <<< "$BOX_Y * 10")")
BOX_Z=$(printf "%.3f" "$(bc -l <<< "$BOX_Z * 10")")
min_ele=$(printf "%.3f" "$(bc -l <<< "$min_ele * 10")")
bulk_xbox=$(printf "%.3f" "$(bc -l <<< "$bulk_xbox * 10")")
bulk_ybox=$(printf "%.3f" "$(bc -l <<< "$bulk_ybox * 10")")
bulk_zbox=$(printf "%.3f" "$(bc -l <<< "$bulk_zbox * 10")")

outfile="${packmol}.inp"
write_packmol_header "$outfile" "2.2" "initial.pdb"

cat >> "$outfile" <<EOF
structure single_electrode.pdb
number 1
fixed  0. 0. 0. 0. 0. 0.
end structure

EOF

echo -e "${OK}Wrote slit Packmol header and electrode block to ${packmol}.inp." | tee -a ./result/b_model.log >&2

mod_max_ele=$(printf "%.3f" "$(bc -l <<< "${max_ele} + 5")")
for i in "${!MOL_arr[@]}"; do
    append_slit_structure_block "$outfile" "${MOL_arr[$i]}" "${MOLnum_arr[$i]}" "$mod_max_ele" "$BOX_Z"
done

echo -e "${OK}Added upper region structures to ${packmol}.inp." | tee -a ./result/b_model.log >&2

mod_min_ele=$(printf "%.3f" "$(bc -l <<< "${min_ele} - 5")")
for i in "${!MOL_arr[@]}"; do
    append_slit_structure_block "$outfile" "${MOL_arr[$i]}" "${MOLnum_arr[$i]}" "1." "$mod_min_ele"
done

echo -e "${OK}Added lower region structures to ${packmol}.inp." | tee -a ./result/b_model.log >&2

outfile="${bulk_pa}.inp"
write_packmol_header "$outfile" "2.0" "initial.pdb"

for i in "${!MOL_arr[@]}"; do
    append_bulk_structure_block "$outfile" "${MOL_arr[$i]}" "${MOLnum_arr_bulk[$i]}"
done

echo -e "${OK}Generated ${bulk_pa}.inp for bulk system." | tee -a ./result/b_model.log >&2
