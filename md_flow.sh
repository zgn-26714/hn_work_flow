#!/bin/bash
rm setting.log


SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
INPUT_FILE="INPUT"

export BLUE='\033[0;34m'
export RED='\033[0;31m'
export GREEN='\033[0;32m'
export YELLOW='\033[1;33m'
export NC='\033[0m'
export OK="${GREEN}[OK]${NC}"
export ERROR="${RED}[ERROR]"

while IFS= read -r line || [[ -n "$line" ]]; do
    line="${line%%#*}"
    line=$(echo "$line" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
    [[ -z "$line" || "$line" =~ ^# ]] && continue
    if [[ "$line" == *"="* ]]; then
        lhs="${line%%=*}"
        rhs="${line#*=}"
        lhs=$(echo "$lhs" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
        rhs=$(echo "$rhs" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
        export "$lhs=$rhs"
        echo "Exported: $lhs=$rhs" >> setting.log
    else
        echo -e "${YELLOW}Warning: Invalid line (no '='): ${line}${NC}" >&2
    fi
done < "$INPUT_FILE"
echo -e "${GREEN}Finish reading INPUT${NC}"

ACTION="${1:-help}"

case "$ACTION" in
    frames)
        echo ">>> Entering get frames branch..."
        setting="${2:-default}"
        case "$setting" in
          default) 
            bash "$SCRIPT_DIR"/build_model/run.sh
            ;;
          bulk)
            export TOP=${bulk_top}
            bash "$SCRIPT_DIR"/build_model/run_for_bulk.sh
          ;;
        esac
        ;;

    run)
        echo ">>> Entering RUN md branch..."
        bash "$SCRIPT_DIR"/run_md/run.sh
        ;;

    analyze)
        echo ">>> Entering analyze branch..."
        bash "$SCRIPT_DIR"/analyze/run.sh $2
        ;;

    # all)
    #     echo ">>> Running ALL steps sequentially..."
    #     # Step 1: Build model
    #     if ! sh build_model.sh; then
    #         echo "❌ Build failed" >&2
    #         exit 1
    #     fi
    #     # Step 2: Equilibrate system
    #     if ! sh equilibration.sh; then
    #         echo "❌ Equilibration failed" >&2
    #         exit 1
    #     fi
    #     # Step 3: Extract frame
    #     if ! sh get_frame.sh; then
    #         echo "❌ Frame extraction failed" >&2
    #         exit 1
    #     fi
    #     # Final success celebration
    #     echo -e "\033[1;32m🎉🎉🎉 ALL STEPS COMPLETED SUCCESSFULLY! 🎉🎉🎉\033[0m"
    #     ;;

    help|*)
        cat <<EOF
================================================================================
                          Molecular Dynamics Pipeline
================================================================================

Description:
  This script orchestrates various stages of a molecular dynamics (MD) workflow.
  It reads configuration from 'INPUT' file and supports modular execution via
  subcommands.

Usage:
  $0 <command> [options...]

Available Commands:

  frames          Generate initial frames for simulation(nvt).
                  → Executes: $SCRIPT_DIR/build_model/run.sh

  run             Submit job to perform MD simulation.
                  → Executes: $SCRIPT_DIR/run_md/run.sh

  analyze         Analyze output data from MD simulations (e.g., energy, force).
                  → Executes: $SCRIPT_DIR/analyze/run.sh

  help            Show this help message (default).

Examples:

  $0 frames
  $0 run
  $0 analyze

Configuration:
  Settings are read from the 'INPUT' file in the current directory.
  Format: KEY = VALUE   (comments with '#' are ignored)

  Example INPUT file:
      TEMPERATURE = 300
      STEPS       = 10000
      OUTPUT_DIR  = ./results

================================================================================
EOF
        exit 0
        ;;
esac


