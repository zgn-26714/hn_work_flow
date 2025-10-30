#!/bin/bash

# set -euo pipefail
rm -rf model
rm -rf bulk
rm -rf build
rm -rf nvt20
rm -rf initial


bash_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

if [[ -f "$GMXRC" ]]; then
    source "$GMXRC" > ./result/b_model.log >&2
else
    echo -e "${ERROR}GMXRC not found: $GMXRC${NC}"  | tee ./result/b_model.log >&2
    exit 1
fi

# 1. 生成单侧电极
if bash ${bash_dir}/source/build_single_electrode.sh; then
    echo -e "${GREEN}"
    echo "##################################################">&2
    echo "  I. success! build_single_electrode.sh completed">&2
    echo "##################################################">&2
    echo -e "${NC}"
else
    echo -e "${ERROR} build_single_electrode.sh failed${NC}" | tee -a ./result/b_model.log  >&2
    exit 1
fi

cp ./model/single_electrode.pdb .
# 2. 调密度
if [ -z "$set_density" ]; then
    output=$(sh "${bash_dir}/get_bulk_density.sh") 

    if [ $? -ne 0 ] || ! echo "$output" | grep -q "^OUTPUT:"; then
        echo -e "${ERROR} get_bulk_density.sh failed${NC}" | tee -a ./result/b_model.log  >&2
        exit 1
    else
        echo -e "${OK}success get bulk density ${output}" | tee -a ./result/b_model.log >&2
        density=$(echo "$output" | grep "^OUTPUT:" | cut -d' ' -f2)
        if ! [[ "$density" =~ ^[+-]?[0-9]*\.?[0-9]+([eE][+-]?[0-9]+)?$ ]]; then
            echo -e "${ERROR}Error: Invalid density value returned: $density${NC}" | tee -a ./result/b_model.log >&2
            echo -e "${ERROR} get_bulk_density.sh failed${NC}" | tee -a ./result/b_model.log  >&2
            exit 1
        else
            echo -e "${OK}success WRITE bulk density ${density}" | tee -a ./result/b_model.log >&2
        fi
    fi
    density=$(printf "%.10f" "$density")
    export set_density="${density}"
fi


if bash ${bash_dir}/build_model.sh; then
    echo -e "${GREEN}"
    echo "##############################################">&2
    echo "  II. success! build_model.sh completed">&2
    echo "  The bulk density has been equilibrated.">&2
    echo "##############################################">&2
    echo -e "${NC}"
else
    echo -e "${ERROR} build_model.sh failed${NC}" | tee -a ./result/b_model.log  >&2
    exit 1
fi

# 3. 平移复制 4. 重新平衡
if bash ${bash_dir}/source/replicate_translate.sh; then
    echo -e "${GREEN}"
    echo "##############################################">&2
    echo "  III. success! replicate_translate.sh completed">&2
    echo "##############################################">&2
    echo -e "${NC}"
else
    echo -e "${ERROR} replicate_translate.sh failed${NC}" | tee -a ./result/b_model.log  >&2
    exit 1
fi


cp ${MDP}.mdp bak_${MDP}.mdp
sed -i -e 's/^freezegrps.*$/freezegrps = CL  CR  GRA/' -e 's/^freezedim.*$/freezedim   = Y Y Y Y Y Y Y Y Y/' ${MDP}.mdp
mv ./build/mini_re.gro ./build/pre_eq.gro

# 5. nvt20
if bash ${bash_dir}/equilibration.sh; then
    echo -e "${GREEN}"
    echo "##############################################">&2
    echo "  IV. success!  equilibration.sh completed">&2
    echo "##############################################">&2
    echo -e "${NC}"
else
    echo -e "${ERROR} equilibration.sh failed${NC}" | tee -a ./result/b_model.log  >&2
    exit 1
fi

# 6. get_frames
if bash ${bash_dir}/get_frames.sh; then
    echo -e "${GREEN}"
    echo "##############################################">&2
    echo "  V. success!  get_frames.sh completed">&2
    echo "##############################################">&2
    echo -e "${NC}"
    rm ./step*.pdb 
    rm ./bulk/*#
    rm ./initial/*#
    rm ./build/*#
    rm ./model/*#
    rm ./nvt20/*#
    echo ""
            echo -e "\033[1;32m🎉🎉🎉  ALL STEPS COMPLETED SUCCESSFULLY!  🎉🎉🎉\033[0m">&2
            echo -e "\033[1;36m      ____  ____  ____  ____  ____  ____       \033[0m">&2
            echo -e "\033[1;36m     ||C ||||O ||||N ||||G ||||R ||||A ||     \033[0m">&2
            echo -e "\033[1;36m     ||__||||__||||__||||__||||__||||__||     \033[0m">&2
            echo -e "\033[1;36m     |/__\||/__\||/__\||/__\||/__\||/__\|     \033[0m">&2
            echo -e "\033[1;33m         ||||T ||||U ||||L ||||A ||||T ||     \033[0m">&2
            echo -e "\033[1;33m         ||||__||||__||||__||||__||||__||     \033[0m">&2
            echo -e "\033[1;33m         ||||/__\||/__\||/__\||/__\||/__\|     \033[0m">&2
            echo -e "\033[1;35m              ||||I ||||O ||||N ||||S ||       \033[0m">&2
            echo -e "\033[1;35m              ||||__||||__||||__||||__||       \033[0m">&2
            echo -e "\033[1;35m              ||||/__\||/__\||/__\||/__\|       \033[0m">&2
            echo -e "\033[1;32m🎉🎉🎉  GO CELEBRATE — YOU EARNED IT!  🎉🎉🎉\033[0m">&2
    echo ""
else
    echo -e "${ERROR} get_frames.sh failed${NC}" | tee -a ./result/b_model.log  >&2
    exit 1
fi
