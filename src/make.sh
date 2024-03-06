#!/bin/bash
GREEN='\033[1;32m'
RED='\033[1;31m'
NC='\033[0m'
if [[ "$#" -eq 1 ]]
then
  err=0
  echo -e "${GREEN}Making ParticleTrackingOF executable for model_$1...${NC}"
	./set_model.sh $1
	make ParticleTrackingOF
  if [[ $? -eq 0 ]]
  then
    mv "../bin/ParticleTrackingOF" "../bin/ParticleTrackingOF_$1"
    echo -e "${GREEN}Done.${NC}"
  else
    echo -e "${RED}Failed.${NC}"
    ((++err))
  fi
  echo -e "${GREEN}Making ParticleTrackingOF_TwoPhaseNonStationary executable for model_$1...${NC}"
  ./set_model.sh $1
  make ParticleTrackingOF_TwoPhaseNonStationary
  if [[ $? -eq 0 ]]
  then
    mv "../bin/ParticleTrackingOF_TwoPhaseNonStationary" "../bin/ParticleTrackingOF_TwoPhaseNonStationary_$1"
    echo -e "${GREEN}Done.${NC}"
  else
    echo -e "${RED}Failed.${NC}"
    ((++err))
  fi
	exit "$err"
fi
echo -e "${RED}Bad parameters.${NC}"
exit 3
