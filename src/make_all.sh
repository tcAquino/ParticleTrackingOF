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

SCRIPT_DIR=$(dirname -- "$( readlink -f -- "$0"; )")
readarray -t models < "$SCRIPT_DIR/models"

err=0
for model in "${models[@]}"
do
    ./make.sh "${1}" "$model" "${2}" ${3}
    ((err+=$?))
done

if [[ err -ne 0 ]] ; then
  echo -e "${RED}$err model compilations failed.${NC}"
  exit 1
else
  echo -e "${GREEN}All models compiled successfuly.${NC}"
  exit 0
fi
