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

if [[ "$#" -ne 2 ]] && [[ "$#" -ne 3 ]] ; then
	echo -e "${RED}Bad parameters.${NC}"
	exit 1
fi

err=0
./make.sh "${1}" advection_diffusion_2d "${2}" ${3}
((err+=$?))
./make.sh "${1}" advection_diffusion_fpt_2d "${2}" ${3}
((err+=$?))
./make.sh "${1}" advection_diffusion_surface_decay_2d "${2}" ${3} 
((err+=$?))
./make.sh "${1}" advection_2d "${2}" ${3}
((err+=$?))
./make.sh "${1}" periodic_cartesian_advection_diffusion_2d "${2}" ${3}
((err+=$?))
./make.sh "${1}" periodic_cartesian_advection_diffusion_fpt_2d "${2}" ${3}
((err+=$?))
./make.sh "${1}" periodic_cartesian_advection_diffusion_surface_decay_2d "${2}" ${3}
((err+=$?))
./make.sh "${1}" periodic_cartesian_advection_2d "${2}" ${3} 
((err+=$?))
./make.sh "${1}" advection_diffusion_3d "${2}" ${3}
((err+=$?))
./make.sh "${1}" advection_diffusion_fpt_3d "${2}" ${3}
((err+=$?))
./make.sh "${1}" advection_diffusion_surface_decay_3d "${2}" ${3}
((err+=$?))
./make.sh "${1}" advection_3d "${2}" ${3}
((err+=$?))
./make.sh "${1}" periodic_cartesian_advection_diffusion_3d "${2}" ${3}
((err+=$?))
./make.sh "${1}" periodic_cartesian_advection_diffusion_fpt_3d "${2}" ${3}
((err+=$?))
./make.sh "${1}" periodic_cartesian_advection_diffusion_surface_decay_3d "${2}" ${3}
((err+=$?))
./make.sh "${1}" periodic_cartesian_advection_3d "${2}" ${3}
((err+=$?))
./make.sh "${1}" bcc_cartesian_advection_diffusion "${2}" ${3}
((err+=$?))
./make.sh "${1}" bcc_cartesian_advection_diffusion_fpt "${2}" ${3}
((err+=$?))
./make.sh "${1}" bcc_cartesian_advection_diffusion_surface_decay "${2}" ${3}
((err+=$?))
./make.sh "${1}" bcc_cartesian_advection "${2}" ${3}
((err+=$?))
./make.sh "${1}" bcc_symmetryplanes_advection_diffusion "${2}" ${3}
((err+=$?))
./make.sh "${1}" bcc_symmetryplanes_advection_diffusion_fpt "${2}" ${3}
((err+=$?))
./make.sh "${1}" bcc_symmetryplanes_advection_diffusion_surface_decay "${2}" ${3}
((err+=$?))
./make.sh "${1}" bcc_symmetryplanes_advection "${2}" ${3}
((err+=$?))

if [[ err -ne 0 ]] ; then
  echo -e "${RED}$err model compilations failed.${NC}"
  exit 1
else
  echo -e "${GREEN}All models compiled successfuly.${NC}"
  exit 0
fi
