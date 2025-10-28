#!/bin/bash

set -euo pipefail

bash_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

if [[ -f "$GMXRC" ]]; then
    source "$GMXRC"
else
    echo -e "${ERROR}GMXRC not found: $GMXRC${NC}"  | tee -a ./result/b_model.log >&2
    exit 1
fi

# 1. 生成单侧电极
if bash ${bash_dir}/source/build_single_electrode.sh; then
    echo -e "${GREEN}"
    echo "##############################################">&2
    echo "  I. success! build_single_electrode.sh completed">&2
    echo "##############################################">&2
    echo -e "${NC}"
else
    echo -e "${ERROR} build_single_electrode.sh failed${NC}" | tee -a ./result/b_model.log  >&2
    exit 1
fi



# 2. 调密度
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

mv build/pre_eq_merge.gro build/pre_eq.gro
# change_top


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
