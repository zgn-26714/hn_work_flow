#!/usr/bin/env python3

import re
import sys
import os
import subprocess
from typing import List


def read_mol_dat(mol_file: str) -> List[str]:
    """
    Read mol.dat file and extract molecule names (one per line).
    Ignore empty lines and comments starting with '#' or ';'.
    """
    if not os.path.exists(mol_file):
        raise FileNotFoundError(f"mol.dat file not found: {mol_file}")

    names = []
    with open(mol_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            stripped = line.strip()
            if not stripped or stripped.startswith('#') or stripped.startswith(';'):
                continue
            parts = stripped.split()
            if not parts:
                continue
            mol_name = parts[0]
            if not mol_name:
                print(f"Warning: Empty molecule name at line {line_num}")
                continue
            names.append(mol_name)
    return names

def modify_packmol(inp_file: str, output_file: str, target_mols: List[str], factor: float):
    """
    Modify the 'number' field in packmol.inp for specified molecules.
    Molecule name is extracted from the .pdb filename in 'structure' block.
    """
    if not os.path.exists(inp_file):
        raise FileNotFoundError(f"Input packmol file not found: {inp_file}")

    with open(inp_file, 'r') as f:
        lines = f.readlines()

    output_lines = []
    in_structure = False
    current_pdb = None

    for line in lines:
        stripped = line.strip()

        # Match start of structure block
        pdb_match = re.match(r'^\s*structure\s+([^#\s]+\.pdb)', stripped, re.IGNORECASE)
        if pdb_match:
            in_structure = True
            current_pdb = pdb_match.group(1)
            output_lines.append(line)
            continue

        # Match end of structure block
        if in_structure and stripped.startswith('end structure'):
            in_structure = False
            current_pdb = None
            output_lines.append(line)
            continue

        # Inside structure block, process 'number'
        if in_structure and stripped.startswith('number'):
            try:
                old_count = int(stripped.split()[1])
            except (IndexError, ValueError) as e:
                print(f"[packmol] Warning: Invalid number line, skipping: {line.strip()}")
                output_lines.append(line)
                continue

            # Extract molecule name from pdb filename (e.g., A.pdb -> A)
            mol_name = re.sub(r'\.pdb$', '', current_pdb, flags=re.IGNORECASE)

            if mol_name in target_mols:
                new_count = int(round(old_count * factor))
                if new_count < 0:
                    raise ValueError(f"Invalid count after scaling: {new_count} for {mol_name}")
                indent = line[:line.find('number')]
                new_line = f"{indent}number {new_count}\n"
                output_lines.append(new_line)
                print(f"[packmol] {mol_name}: {old_count} -> {new_count} (x{factor:.3f})")
            else:
                output_lines.append(line)
        else:
            output_lines.append(line)

    with open(output_file, 'w') as f:
        f.writelines(output_lines)
    print(f"Modified packmol input written to: {output_file}")


def modify_topology(top_file: str, output_file: str, target_mols: List[str], factor: float):
    """
    Modify the molecule counts in [ molecules ] section of the topology (.top) file.
    """
    if not os.path.exists(top_file):
        raise FileNotFoundError(f"Topology file not found: {top_file}")

    with open(top_file, 'r') as f:
        lines = f.readlines()

    output_lines = []
    in_molecules = False

    for line in lines:
        stripped = line.strip()

        if stripped.startswith('[ molecules ]'):
            in_molecules = True
            output_lines.append(line)
            continue

        # End of [ molecules ] section when another directive starts
        if in_molecules and stripped.startswith('[') and 'molecules' not in stripped.lower():
            in_molecules = False

        if in_molecules:
            if not stripped or stripped.startswith(';'):
                output_lines.append(line)
                continue

            parts = stripped.split()
            if len(parts) >= 2:
                mol_name = parts[0]
                try:
                    old_count = int(parts[1])
                except ValueError:
                    print(f"[top] Warning: Invalid count in line: {line.strip()}")
                    output_lines.append(line)
                    continue

                if mol_name in target_mols:
                    new_count = int(round(old_count * factor))
                    if new_count < 0:
                        raise ValueError(f"Invalid count after scaling: {new_count} for {mol_name}")
                    # Preserve formatting as much as possible
                    new_line = f"  {mol_name:8s} {new_count:>8d}\n"
                    output_lines.append(new_line)
                    print(f"[top] {mol_name}: {old_count} -> {new_count} (x{factor:.3f})")
                else:
                    output_lines.append(line)
            else:
                output_lines.append(line)
        else:
            output_lines.append(line)

    with open(output_file, 'w') as f:
        f.writelines(output_lines)
    print(f"Modified topology written to: {output_file}")


def main():
    script_name = os.path.basename(__file__)  # 只取文件名，如 "script.py"
    iter = os.environ['iter']
    # 或者用 full path: __file__

    print(f"""
===================================================================================================
 ---------------------------------------- iter="{iter}" -----------------------------------------
【INFO】Starting pre-equilibration and density analysis workflow
Script    : {script_name}
===================================================================================================""")

    inp_file = os.environ['INP_FILE']
    top_file = os.environ['TOP']
    mol_name = os.environ['MOL_name']
    target_molecules = mol_name.split()
    # Read scaling factor from environment variable

    try:
        factor = float(sys.argv[1])
        if factor < 0:
            print("Error: SCALE_FACTOR must be non-negative")
            exit(1)
    except (IndexError, ValueError):
        print("Error: Please provide a valid non-negative scaling factor as the first argument")
        exit(1)

    # Modify packmol input
    try:
        modify_packmol(inp_file, inp_file, target_molecules, factor)
    except Exception as e:
        print(f"Error: Failed to modify packmol.inp: {e}")
        exit(1)

    # Modify topology file
    top_file = f"{top_file}.top"
    try:
        modify_topology(top_file, top_file, target_molecules, factor)
    except Exception as e:
        print(f"Error: Failed to modify topology file: {e}")
        exit(1)

    print("✅ All files modified successfully!")


if __name__ == '__main__':
    main()