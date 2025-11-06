#!/bin/bash
set -euo pipefail  

# 
make_job() {
    local start_case=$1
    local end_case=$2
    local rundir=$3
    local workdir=$4
    if [[ "${isRerun}" -eq 1 ]]; then
        local objjob="$workdir/jobfile/rerun_j${V}V${ic}pscase${start_case}-${end_case}.job"
    else
        local objjob="$workdir/jobfile/j${V}V${ic}pscase${start_case}-${end_case}.job"
    fi
    mkdir -p "$workdir/jobfile"
    local SCRIPT_PATH
    SCRIPT_PATH="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
    
    case "$server_machine" in
       "eninstein") cp "$SCRIPT_PATH/eninstein.job" "$objjob" ;;
        "4090")     cp "$SCRIPT_PATH/4090.job" "$objjob" ;; 
        "wuchao")   cp "$SCRIPT_PATH/wuchao.job" "$objjob" ;;
    esac
    # 
    
    case "$server_machine" in
    "eninstein") 
        sed -i "
            2s|.*|#PBS -q $queue|;
            3s|.*|#PBS -N chargingcase${start_case}-${end_case}${V}V${ic}|;
            /k=/s|.*|k=$start_case|;
            /h=/s|.*|h=$num|;
            /workdir=/s|.*|workdir=$rundir|;
            /mode=/s|.*|mode=$mode|;
            /isRerun=/s|.*|isRerun=$isRerun|;
            /name=/s|.*|name=$DEFFNM|
        " "$objjob" 
        sed -i 's|.*export ONFLY_DENSITY3D=.*|export ONFLY_DENSITY3D="-low ['"${LOW_ONFLY}"'] -up ['"${UP_ONFLY}"'] -nbin ['"${NBIN_ONFLY}"'] -n index.ndx -sel ['"${MOL_name}"'] -calc '"${MODE_ONFLY}"'"|' "$objjob"
        case "$queue" in
            short)  sed -i '6s/.*/#PBS -l walltime=36:00:00/' "$objjob" ;;
            new)    sed -i '5s/.*/#PBS -l nodes=1:ppn=64/' "$objjob" ;;
            fast)   sed -i '5s/.*/#PBS -l nodes=1:ppn=16/' "$objjob" ;;
            long)   sed -i '6s/.*/#PBS -l walltime=720:00:00/' "$objjob" ;;
        esac

        if [[ "${isRerun}" -eq 1 ]]; then
            cd "$rundir"/case"$start_case"/rerun_case
            sed -i "/casedir=/s|.*|casedir=\$workdir/case\$k/rerun_case|" "$objjob"
        else
            cd "$rundir"/case"$start_case"
            pwd
        fi
        qsub "$objjob"
        cd -
        ;;
    "4090") 
        sed -i "
            5s|.*|#SBATCH --partition=${queue}|;
            2s|.*|#SBATCH --job-name=chargingcase${start_case}-${end_case}${V}V${ic}|;
            /k=/s|.*|k=$start_case|;
            /h=/s|.*|h=$num|;
            /workdir=/s|.*|workdir=$rundir|;
            /mode=/s|.*|mode=$mode|;
            /isRerun=/s|.*|isRerun=$isRerun|;
            /name=/s|.*|name=$DEFFNM|
        " "$objjob" 
        sed -i 's|.*export ONFLY_DENSITY3D=.*|export ONFLY_DENSITY3D="-low ['"${LOW_ONFLY}"'] -up ['"${UP_ONFLY}"'] -nbin ['"${NBIN_ONFLY}"'] -n index.ndx -sel ['"${MOL_name}"'] -calc '"${MODE_ONFLY}"'"|' "$objjob"
        if [[ "${isRerun}" -eq 1 ]]; then
            cd "$rundir"/case"$start_case"/rerun_case
            sed -i "/casedir=/s|.*|casedir=\$workdir/case\$k/rerun_case|" "$objjob"
        else
            cd "$rundir"/case"$start_case"
        fi
        sbatch "$objjob"
        cd -
        ;;
    "wuchao") 
        sed -i "
            7s|.*|#DSUB -n chargingcase${start_case}-${end_case}${V}V${ic}|;
            /k=/s|.*|k=$start_case|;
            /h=/s|.*|h=$num|;
            /workdir=/s|.*|workdir=$rundir|;
            /mode=/s|.*|mode=$mode|;
            /isRerun=/s|.*|isRerun=$isRerun|;
            /name=/s|.*|name=$DEFFNM|
        " "$objjob" 
        sed -i 's|.*export ONFLY_DENSITY3D=.*|export ONFLY_DENSITY3D="-low ['"${LOW_ONFLY}"'] -up ['"${UP_ONFLY}"'] -nbin ['"${NBIN_ONFLY}"'] -n index.ndx -sel ['"${MOL_name}"'] -calc '"${MODE_ONFLY}"'"|' "$objjob"
        chmod +x $objjob
        if [[ "${isRerun}" -eq 1 ]]; then
            cd "$rundir"/case"$start_case"/rerun_case
            sed -i "/casedir=/s|.*|casedir=\$workdir/case\$k/rerun_case|" "$objjob"
        else
            cd "$rundir"/case"$start_case"
        fi
        dsub -s $objjob
        cd -
        ;;
    *)  
        echo "❌ Unknown server machine: $server_machine" | tee -a ./result/run_md.log >&2
        return 1
        ;;
esac
}

# 
main() {
    # 
    workdir="$(pwd)"
    echo $workdir |tee  /data1/huangnan/PC/charging/298k/2V/0ps/debug 
    local current_cas=${START}
    local batch_num=1
    
    while (( current_cas <= ${END} )); do
        local batch_end=$((current_cas + num - 1))
        if (( batch_end > ${END} )); then
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
        ((batch_num++))
    done


    # Summary
    echo -e "${GREEN}✅ Case generation completed successfully!${NC}" | tee -a ./result/run_md.log >&2
}

main "$@"