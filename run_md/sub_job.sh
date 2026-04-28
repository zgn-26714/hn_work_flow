#!/bin/bash
set -euo pipefail

SCRIPT_PATH="$(builtin cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"

sed_in_place() {
    if sed --version >/dev/null 2>&1; then
        sed -i "$@"
    else
        sed -i '' "$@"
    fi
}

job_file_path() {
    local start_case="$1"
    local end_case="$2"
    local workdir="$3"

    if [[ "${isRerun}" -eq 1 ]]; then
        printf '%s/jobfile/rerun_j%sV%spscase%s-%s.job\n' \
            "$workdir" "$V" "$ic" "$start_case" "$end_case"
    else
        printf '%s/jobfile/j%sV%spscase%s-%s.job\n' \
            "$workdir" "$V" "$ic" "$start_case" "$end_case"
    fi
}

job_template_path() {
    case "$server_machine" in
        eninstein) printf '%s/eninstein.job\n' "$SCRIPT_PATH" ;;
        4090)      printf '%s/4090.job\n' "$SCRIPT_PATH" ;;
        wuchao)    printf '%s/wuchao.job\n' "$SCRIPT_PATH" ;;
        jiaocha)   printf '%s/jiaocha.job\n' "$SCRIPT_PATH" ;;
        *)
            echo "❌ Unknown server machine: $server_machine" | tee -a ./result/run_md.log >&2
            return 1
            ;;
    esac
}

build_onfly_density3d() {
    printf 'export ONFLY_DENSITY3D="-low [%s] -up [%s] -nbin [%s] -n index.ndx -sel [%s] -calc %s"' \
        "$LOW_ONFLY" "$UP_ONFLY" "$NBIN_ONFLY" "$MOL_name_run" "$MODE_ONFLY"
}

apply_common_job_settings() {
    local objjob="$1"
    local start_case="$2"
    local rundir="$3"
    local onfly_density3d

    onfly_density3d="$(build_onfly_density3d)"
    sed_in_place \
        -e "/^[[:space:]]*k[[:space:]]*=/s|.*|k=$start_case|" \
        -e "/^[[:space:]]*h[[:space:]]*=/s|.*|h=$num|" \
        -e "/^[[:space:]]*workdir[[:space:]]*=/s|.*|workdir=$rundir|" \
        -e "/^[[:space:]]*mode[[:space:]]*=/s|.*|mode=$mode|" \
        -e "/^[[:space:]]*isRerun[[:space:]]*=/s|.*|isRerun=$isRerun|" \
        -e "/^[[:space:]]*name[[:space:]]*=/s|.*|name=$DEFFNM|" \
        -e "/^[[:space:]]*mdbin[[:space:]]*=/s|.*|mdbin=$GMX_run|" \
        -e "s|export ONFLY_FLAGS=0|export ONFLY_FLAGS=$ONFLY_FLAGS|" \
        -e "s|.*export ONFLY_DENSITY3D=.*|$onfly_density3d|" \
        "$objjob"
}

configure_eninstein_job() {
    local objjob="$1"
    local start_case="$2"
    local end_case="$3"
    local rundir="$4"
    local walltime

    apply_common_job_settings "$objjob" "$start_case" "$rundir"
    sed_in_place \
        -e "2s|.*|#PBS -q $queue|" \
        -e "3s|.*|#PBS -N chargingcase${start_case}-${end_case}${V}V${ic}|" \
        -e "5s|.*|#PBS -l nodes=1:ppn=${server_core}|" \
        "$objjob"

    case "$queue" in
        short) walltime="36:00:00" ;;
        long)  walltime="720:00:00" ;;
        *)     walltime="" ;;
    esac

    if [[ -n "$walltime" ]]; then
        sed_in_place -e "6s|.*|#PBS -l walltime=${walltime}|" "$objjob"
    fi
}

configure_4090_job() {
    local objjob="$1"
    local start_case="$2"
    local end_case="$3"
    local rundir="$4"

    apply_common_job_settings "$objjob" "$start_case" "$rundir"
    sed_in_place \
        -e "2s|.*|#SBATCH --job-name=chargingcase${start_case}-${end_case}${V}V${ic}|" \
        -e "5s|.*|#SBATCH --partition=${queue}|" \
        -e "s|^[[:space:]]*nofp[[:space:]]*=[[:space:]]*[0-9][0-9]*;[[:space:]]*$|nofp=${server_core};|" \
        "$objjob"
}

configure_wuchao_job() {
    local objjob="$1"
    local start_case="$2"
    local end_case="$3"
    local rundir="$4"

    apply_common_job_settings "$objjob" "$start_case" "$rundir"
    sed_in_place \
        -e "7s|.*|#DSUB -n chargingcase${start_case}-${end_case}${V}V${ic}|" \
        -e "s|^#DSUB -R cpu=[0-9][0-9]*;mem=|#DSUB -R cpu=${server_core};mem=|" \
        -e "s|^[[:space:]]*nofp[[:space:]]*=[[:space:]]*[0-9][0-9]*;[[:space:]]*$|nofp=${server_core};|" \
        "$objjob"
}

configure_jiaocha_job() {
    local objjob="$1"
    local start_case="$2"
    local end_case="$3"
    local rundir="$4"

    apply_common_job_settings "$objjob" "$start_case" "$rundir"
    sed_in_place \
        -e "2s|.*|#SBATCH --job-name=chargingcase${start_case}-${end_case}${V}V${ic}|" \
        -e "5s|.*|#SBATCH --ntasks-per-node=${server_core} ### 每个节点所运行的进程数为${server_core}|" \
        -e "s|-ntomp 56|-ntomp ${server_core}|g" \
        "$objjob"
}

configure_job_template() {
    local objjob="$1"
    local start_case="$2"
    local end_case="$3"
    local rundir="$4"

    case "$server_machine" in
        eninstein) configure_eninstein_job "$objjob" "$start_case" "$end_case" "$rundir" ;;
        4090)      configure_4090_job "$objjob" "$start_case" "$end_case" "$rundir" ;;
        wuchao)    configure_wuchao_job "$objjob" "$start_case" "$end_case" "$rundir" ;;
        jiaocha)   configure_jiaocha_job "$objjob" "$start_case" "$end_case" "$rundir" ;;
        *)
            echo "❌ Unknown server machine: $server_machine" | tee -a ./result/run_md.log >&2
            return 1
            ;;
    esac
}

submit_case_dir() {
    local rundir="$1"
    local start_case="$2"

    if [[ "${isRerun}" -eq 1 ]]; then
        printf '%s/case%s/rerun_case\n' "$rundir" "$start_case"
    else
        printf '%s/case%s\n' "$rundir" "$start_case"
    fi
}

apply_rerun_case_dir() {
    local objjob="$1"

    if [[ "${isRerun}" -eq 1 ]]; then
        sed_in_place -e '/casedir=/s|.*|casedir=\$workdir/case\$k/rerun_case|' "$objjob"
    fi
}

submit_job_file() {
    local objjob="$1"
    local submit_dir="$2"
    local previous_dir submit_rc

    previous_dir="$(pwd)"
    cd "$submit_dir"

    case "$server_machine" in
        eninstein)
            set +e
            qsub "$objjob"
            submit_rc=$?
            set -e
            ;;
        4090)
            set +e
            sbatch "$objjob"
            submit_rc=$?
            set -e
            ;;
        wuchao)
            chmod +x "$objjob"
            set +e
            dsub -s "$objjob"
            submit_rc=$?
            set -e
            ;;
        jiaocha)
            chmod +x "$objjob"
            set +e
            sbatch -s "$objjob"
            submit_rc=$?
            set -e
            ;;
        *)
            cd "$previous_dir"
            echo "❌ Unknown server machine: $server_machine" | tee -a ./result/run_md.log >&2
            return 1
            ;;
    esac
    cd "$previous_dir"
    return "$submit_rc"
}

make_job() {
    local start_case="$1"
    local end_case="$2"
    local rundir="$3"
    local workdir="$4"
    local objjob template_path submit_dir

    objjob="$(job_file_path "$start_case" "$end_case" "$workdir")"
    template_path="$(job_template_path)"
    submit_dir="$(submit_case_dir "$rundir" "$start_case")"

    mkdir -p "$workdir/jobfile"
    cp "$template_path" "$objjob"

    configure_job_template "$objjob" "$start_case" "$end_case" "$rundir"
    apply_rerun_case_dir "$objjob"
    submit_job_file "$objjob" "$submit_dir"
}

main() {
    local workdir current_cas batch_end rundir

    workdir="$(pwd)"
    current_cas=${START}

    while (( current_cas <= END )); do
        batch_end=$((current_cas + num - 1))
        if (( batch_end > END )); then
            batch_end=${END}
        fi

        if [[ "${isRerun}" -eq 1 ]]; then
            rundir="$workdir"
        else
            rundir="$workdir/charging/${Temperature}k/${V}V/${ic}ps"
        fi

        echo "📦 Creating job for cases ${current_cas}-${batch_end}..." | tee -a ./result/run_md.log >&2
        make_job "$current_cas" "$batch_end" "$rundir" "$workdir"
        current_cas=$((batch_end + 1))
    done

    echo -e "${GREEN:-}✅ Case generation completed successfully!${NC:-}" | tee -a ./result/run_md.log >&2
}

main "$@"
