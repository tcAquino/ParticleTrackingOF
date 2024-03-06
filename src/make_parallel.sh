#!/bin/bash
GREEN='\033[1;32m'
RED='\033[1;31m'
NC='\033[0m'
if [[ "$#" -eq 1 ]]
then
  err=0
  echo -e "${GREEN}Making ParticleTrackingOF executable for model_$1_parallel...${NC}"
	./set_model_parallel.sh $1
	make ParticleTrackingOF_Parallel
  if [[ $? -eq 0 ]]
  then
    mv "../bin/ParticleTrackingOF_Parallel" "../bin/ParticleTrackingOF_$1_parallel${NC}"
    echo -e "${GREEN}Done.${NC}"
  else
    echo -e "${RED}Failed.${NC}"
    ((++err))
  fi
  echo -e "${GREEN}Making ParticleTrackingOF_TwoPhaseNonStationary executable for model_$1_parallel...${NC}"
  ./set_model_parallel.sh $1
  make ParticleTrackingOF_TwoPhaseNonStationary_Parallel
  if [[ $? -eq 0 ]]
  then
    mv "../bin/ParticleTrackingOF_TwoPhaseNonStationary_Parallel" "../bin/ParticleTrackingOF_TwoPhaseNonStationary_$1_parallel"
    echo -e "${GREEN}Done.${NC}"
  else
    echo -e "${RED}Failed.${NC}"
    ((++err))
  fi
	exit "$err"
fi
echo "${RED}Bad parameters.${NC}"
exit 2
