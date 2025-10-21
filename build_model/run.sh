#!/bin/bash

#  build_model.sh
rm -rf bulk
rm -rf build
rm -rf nvt20
rm -rf initial

bash_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
rm -rf result
mkdir -p result

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



if sh "${bash_dir}/build_model.sh"; then
    echo -e "${GREEN}"
    echo "##############################################">&2
    echo "  I. success! build_model.sh completed">&2
    echo "  The bulk density has been equilibrated.">&2
    echo "##############################################">&2
    echo -e "${NC}"
    rm step*.pdb
    # 执行 equilibration.sh
    if sh "${bash_dir}/equilibration.sh"; then
        echo -e "${GREEN}"
        echo "##############################################">&2
        echo "  II. success! equilibration.sh completed">&2
        echo "##############################################">&2
        echo -e "${NC}"
        # 执行 get_frame.sh
        if sh "${bash_dir}/get_frames.sh"; then
            cp "${TOP}".top result.top
            mv baktop1.top "${TOP}".top
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
            echo -e "${ERROR} get_frame.sh failed${NC}"  | tee -a ./result/b_model.log >&2
            exit 1
        fi

    else
        echo -e "${ERROR} equilibration.sh failed${NC}"  | tee -a ./result/b_model.log >&2
        exit 1
    fi

else
    echo -e "${ERROR} build_model.sh failed${NC}" | tee -a ./result/b_model.log  >&2
    exit 1
fi