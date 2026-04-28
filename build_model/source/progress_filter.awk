#!/usr/bin/awk -f
# Filter GROMACS mdrun progress output for real-time in-place display.
# Splits records on \r so each progress update is processed immediately.

BEGIN {
    RS = "\r|\n"
}

/^step/ {
    # in-place update: clear line, print, return to start
    printf "\r\033[K%s", $0 > "/dev/stderr"
    fflush("/dev/stderr")
}

/^Performance/ {
    # final summary on its own line
    printf "\n%s\n", $0 > "/dev/stderr"
    fflush("/dev/stderr")
}
