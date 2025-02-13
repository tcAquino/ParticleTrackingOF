#!/bin/bash
# Parameters:
# - Main cpp file name, without extension
# - Model name
# - Parallelization option (serial or parallel)
# - Build mode (release or debug, release by default)
# Exit codes:
# 0 - Success
# 1 - Failure

GREEN='\033[1;32m'
RED='\033[1;31m'
NC='\033[0m'
BUILD=release
PTHBIN_BASE="../bin"

if [[ "$#" -eq 4 ]] ; then
	if [[ "$4" != "debug" ]] && [[ "$4" != "release" ]] ; then
		echo -e "${RED}Build mode should be release or debug, got $4.${NC}"
		exit 1
	fi
	BUILD="$4"
fi

if [[ "$#" -ne 3 ]] && [[ "$#" -ne 4 ]] ; then
	echo -e "${RED}Bad parameters.${NC}"
	exit 1
fi

./set_model.sh "${1}" "${2}" "${3}"  
model_set=$?
if [[ "${model_set}" -ne 0 ]] ; then
	exit 1
fi

file_modifier="";
if [[ ${3,,} = "parallel" ]] ; then
	file_modifier="_parallel"
fi

echo -e "${GREEN}Making ${1} ${3,,} executable for model ${2} ($BUILD build)...${NC}"
make "BUILD=$BUILD" "PTHBIN_BASE=$PTHBIN_BASE" "${1}${file_modifier}"
if [[ $? -eq 0 ]] ; then
	mv "$PTHBIN_BASE/$BUILD/${1}${file_modifier}" "$PTHBIN_BASE/$BUILD/${1}${file_modifier}_${2}"
	if [[ "$?" -ne 0 ]] ; then
		echo -e "${RED}Failed (compilation suceeded, but move failed).${NC}"
		exit 1
	fi
	echo -e "${GREEN}Done.${NC}"
else
	echo -e "${RED}Failed.${NC}"
	exit 1
fi

exit 0
