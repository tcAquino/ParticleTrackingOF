#!/bin/bash
# Parameters:
# - Main cpp file name, without extension
# - Parallelization option
# - Build mode (optional)
# Exit codes:
# 0 - Success
# 1 - Failure

GREEN='\033[1;32m'
RED='\033[1;31m'
NC='\033[0m'

SCRIPT_DIR=$(dirname -- "$( readlink -f -- "$0"; )")
readarray -t a < "SCRIPT_DIR/models"

declare -a models=("advection_diffusion_2d"
		   "advection_diffusion_fpt_2d",
		   "advection_diffusion_surface_decay_2d"
		   "advection_diffusion_surface_order2_2d"
		   "advection_2d"
		   "periodic_cartesian_advection_diffusion_2d"
		   "periodic_cartesian_advection_diffusion_fpt_2d"
		   "periodic_cartesian_advection_diffusion_surface_decay_2d"
		   "periodic_cartesian_advection_diffusion_surface_order2_2d"
		   "periodic_cartesian_advection_2d"
		   "advection_diffusion_3d"
		   "advection_diffusion_fpt_3d"
		   "advection_diffusion_surface_decay_3d"
		   "advection_diffusion_surface_order2_3d"
		   "advection_3d"
		   "periodic_cartesian_advection_diffusion_3d"
		   "periodic_cartesian_advection_diffusion_fpt_3d"
		   "periodic_cartesian_advection_diffusion_surface_decay_3d"
		   "periodic_cartesian_advection_diffusion_surface_order2_3d"
		   "periodic_cartesian_advection_3d"
		   "bcc_cartesian_advection_diffusion"
		   "bcc_cartesian_advection_diffusion_fpt"
		   "bcc_cartesian_advection_diffusion_surface_decay"
		   "bcc_cartesian_advection_diffusion_surface_order2"
		   "bcc_cartesian_advection"
		   "bcc_symmetryplanes_advection_diffusion"
		   "bcc_symmetryplanes_advection_diffusion_fpt"
		   "bcc_symmetryplanes_advection_diffusion_surface_decay"
		   "bcc_symmetryplanes_advection_diffusion_surface_order2"
		   "bcc_symmetryplanes_advection"
		  )

if [[ "$#" -ne 2 ]] && [[ "$#" -ne 3 ]] ; then
    echo -e "${RED}Bad parameters.${NC}"
    exit 1
fi

err=0
for i in "${models[@]}"
do
    ./make.sh "${1}" advection_diffusion_2d "${2}" ${3}
    ((err+=$?))
done

if [[ err -ne 0 ]] ; then
  echo -e "${RED}$err model compilations failed.${NC}"
  exit 1
else
  echo -e "${GREEN}All models compiled successfuly.${NC}"
  exit 0
fi
