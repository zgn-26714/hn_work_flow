# bash completion for md_flow.sh
# Usage: source this file in your ~/.bashrc
#   echo "source /path/to/bash_completion.sh" >> ~/.bashrc

_md_flow() {
    local cur prev words cword
    _init_completion || return

    # top-level subcommands
    local commands="frames run analyze d_data help"

    if [[ $cword -eq 1 ]]; then
        COMPREPLY=($(compgen -W "$commands" -- "$cur"))
        return
    fi

    case "${words[1]}" in
        frames)
            if [[ $cword -eq 2 ]]; then
                COMPREPLY=($(compgen -W "default bulk slit clear" -- "$cur"))
            elif [[ "${words[2]:-}" == "clear" ]]; then
                COMPREPLY=($(compgen -W "--dry-run --keep-bin --help" -- "$cur"))
            fi
            ;;
        run)
            if [[ $cword -eq 2 ]]; then
                COMPREPLY=($(compgen -W "default rerun" -- "$cur"))
            fi
            ;;
        analyze)
            if [[ $cword -eq 2 ]]; then
                COMPREPLY=($(compgen -W "eleQ onfly mdheat MD_PP" -- "$cur"))
            fi
            ;;
        d_data)
            if [[ $cword -eq 2 ]]; then
                COMPREPLY=($(compgen -W "default" -- "$cur"))
            fi
            ;;
    esac
}

complete -F _md_flow md_flow.sh
complete -F _md_flow ./md_flow.sh
