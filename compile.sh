#!/bin/bash

#./src/conda/conda_env_setup.sh --all --conda_root /home/zyin/miniconda3

##new compilation
sed -i 's/<USER>/'$USER'/' src/default_LoCo.jsonc
#./src/conda/conda_env_setup.sh -e base --conda_root /lustre/f2/dev/Zun.Yin/anaconda --env_dir /lustre/f2/dev/Zun.Yin/anaconda/envs
#./src/conda/conda_env_setup.sh -e python3_base --conda_root /lustre/f2/dev/Zun.Yin/anaconda --env_dir /lustre/f2/dev/Zun.Yin/anaconda/envs
./src/conda/conda_env_setup.sh -e R_base --env_dir /nbhome/$USER/.conda/envs
# It is still unclear to me why the generated mdtf file includes the wrong envs path, this fixes it
sed -i 's|envs/_MDTF_base|envs/_MDTF_R_base|' mdtf
