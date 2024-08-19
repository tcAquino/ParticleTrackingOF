#!/bin/bash
# Exit code: number of models that failed to compile

GREEN='\033[1;32m'
RED='\033[1;31m'
NC='\033[0m'
BUILD=release

if [[ "$#" -eq 2 ]] ; then
	if [[ "$2" != "debug" ]] && [[ "$2" != "release" ]] ; then
		echo -e "${RED}Build mode should be release or debug, got $2.${NC}"
	       	exit 3
	fi
	BUILD="$2"
fi

PTHBIN_BASE="../bin"
if [[ "$#" -eq 1 ]] || [[ "$#" -eq 2 ]] ; then
	err=0

	./set_model_parallel.sh $1
	model_set=$?

      	echo -e "${GREEN}Making ParticleTrackingOF executable for model_$1_parallel ($BUILD build)...${NC}"
	[[ "$model_set" -ne 1 ]] && [[ "$model_set" -lt 3 ]] && make "BUILD=$BUILD" "PTHBIN_BASE=$PTHBIN_BASE" ParticleTrackingOF_Parallel
	if [[ $? -eq 0 ]] ; then
		mv "$PTHBIN_BASE/$BUILD/ParticleTrackingOF_Parallel" "$PTHBIN_BASE/$BUILD/ParticleTrackingOF_$1_parallel"
		echo -e "${GREEN}Done.${NC}"
	else
		echo -e "${RED}Failed.${NC}"
		((++err))
	fi
	
	echo -e "${GREEN}Making ParticleTrackingOF_TwoPhaseNonStationary executable for model_$1_parallel ($BUILD build)...${NC}"
	[[ "$model_set" -ne 2 ]] && [[ "$model_set" -lt 3 ]] && make "BUILD=$BUILD" "PTHBIN_BASE=$PTHBIN_BASE" ParticleTrackingOF_TwoPhaseNonStationary_Parallel
	if [[ $? -eq 0 ]] ; then
		mv "$PTHBIN_BASE/$BUILD/ParticleTrackingOF_TwoPhaseNonStationary_Parallel" "$PTHBIN_BASE/$BUILD/ParticleTrackingOF_TwoPhaseNonStationary_$1_parallel"
		echo -e "${GREEN}Done.${NC}"
	else
		echo -e "${RED}Failed.${NC}"
	       	((++err))
	fi

	exit "$err"
fi

echo -e "${RED}Bad parameters.${NC}"
exit 2
