#!/bin/bash
GREEN='\033[1;32m'
RED='\033[1;31m'
NC='\033[0m'
err=0
BUILD=release
if [[ "$#" -eq 1 ]]
  then
  if [[ "$1" != "debug" ]] && [[ "$1" != "release" ]]
  then
    echo -e "${RED}Build mode should be release or debug, got $1.${NC}"
    exit 3
  fi
  BUILD="$1"
fi
./make.sh advection_diffusion_2d "$BUILD"
((err+=$?))
./make.sh advection_diffusion_fpt_2d "$BUILD"
((err+=$?))
./make.sh advection_diffusion_surface_decay_2d "$BUILD" 
((err+=$?))
./make.sh advection_2d "$BUILD"
((err+=$?))
./make.sh periodic_cartesian_advection_diffusion_2d "$BUILD"
((err+=$?))
./make.sh periodic_cartesian_advection_diffusion_fpt_2d "$BUILD"
((err+=$?))
./make.sh periodic_cartesian_advection_diffusion_surface_decay_2d "$BUILD"
((err+=$?))
./make.sh periodic_cartesian_advection_2d "$BUILD" 
((err+=$?))
./make.sh advection_diffusion_3d "$BUILD"
((err+=$?))
./make.sh advection_diffusion_fpt_3d "$BUILD"
((err+=$?))
./make.sh advection_diffusion_surface_decay_3d "$BUILD"
((err+=$?))
./make.sh advection_3d "$BUILD"
((err+=$?))
./make.sh periodic_cartesian_advection_diffusion_3d "$BUILD"
((err+=$?))
./make.sh periodic_cartesian_advection_diffusion_fpt_3d "$BUILD"
((err+=$?))
./make.sh periodic_cartesian_advection_diffusion_surface_decay_3d "$BUILD"
((err+=$?))
./make.sh periodic_cartesian_advection_3d "$BUILD"
((err+=$?))
./make.sh bcc_cartesian_advection_diffusion "$BUILD"
((err+=$?))
./make.sh bcc_cartesian_advection_diffusion_fpt "$BUILD"
((err+=$?))
./make.sh bcc_cartesian_advection_diffusion_surface_decay "$BUILD"
((err+=$?))
./make.sh bcc_cartesian_advection "$BUILD"
((err+=$?))
./make.sh bcc_symmetryplanes_advection_diffusion "$BUILD"
((err+=$?))
./make.sh bcc_symmetryplanes_advection_diffusion_fpt "$BUILD"
((err+=$?))
./make.sh bcc_symmetryplanes_advection_diffusion_surface_decay "$BUILD"
((err+=$?))
./make.sh bcc_symmetryplanes_advection "$BUILD"
((err+=$?))
if [[ err -ne 0 ]]
then
  echo -e "${RED}$err model compilations failed.${NC}"
  exit 1
else
  echo -e "${GREEN}All models compiled successfuly.${NC}"
  exit 0
fi
