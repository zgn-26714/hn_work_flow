#!/bin/bash

#  build_model.sh

bash_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

if sh "${bash_dir}/build_model.sh"; then
    echo "##############################################"
    echo "  I. success! build_model.sh completed"
    echo "  The bulk density has been equilibrated."
    echo "##############################################"

    # 执行 equilibration.sh
    if sh "${bash_dir}/equilibration.sh"; then
        echo "##############################################"
        echo "  II. success! equilibration.sh completed"
        echo "##############################################"
        # 执行 get_frame.sh
        if sh "${bash_dir}/get_frames.sh"; then
            echo ""
            echo -e "\033[1;32m🎉🎉🎉  ALL STEPS COMPLETED SUCCESSFULLY!  🎉🎉🎉\033[0m"
            echo -e "\033[1;36m      ____  ____  ____  ____  ____  ____       \033[0m"
            echo -e "\033[1;36m     ||C ||||O ||||N ||||G ||||R ||||A ||     \033[0m"
            echo -e "\033[1;36m     ||__||||__||||__||||__||||__||||__||     \033[0m"
            echo -e "\033[1;36m     |/__\||/__\||/__\||/__\||/__\||/__\|     \033[0m"
            echo -e "\033[1;33m         ||||T ||||U ||||L ||||A ||||T ||     \033[0m"
            echo -e "\033[1;33m         ||||__||||__||||__||||__||||__||     \033[0m"
            echo -e "\033[1;33m         ||||/__\||/__\||/__\||/__\||/__\|     \033[0m"
            echo -e "\033[1;35m              ||||I ||||O ||||N ||||S ||       \033[0m"
            echo -e "\033[1;35m              ||||__||||__||||__||||__||       \033[0m"
            echo -e "\033[1;35m              ||||/__\||/__\||/__\||/__\|       \033[0m"
            echo -e "\033[1;32m🎉🎉🎉  GO CELEBRATE — YOU EARNED IT!  🎉🎉🎉\033[0m"
            echo ""
        else
            echo "ERROR: get_frame.sh failed" >&2
            exit 1
        fi

    else
        echo "ERROR: equilibration.sh failed" >&2
        exit 1
    fi

else
    echo "ERROR: build_model.sh failed" >&2
    exit 1
fi