#!/bin/bash
# Parameters:
# - Main cpp file name, without extension
# - Model name
# - Parallelization option (serial or parallel)
# Exit codes:
# 0 - Success
# 1 - Failure

if [[ "$#" -ne 3 ]] ; then
    echo "$(basename "$0") : Expected model name and parallelization option."
    exit 1
fi

if [[ "${3,,}" = "serial" ]] ; then
    parallel="Serial"
elif [[ "${3,,}" = "parallel" ]] ; then
    parallel="Parallel"
else
    echo "$(basename "$0") : Parallelization option should be serial or parallel, got ${3}."
    exit 1
fi

match="$(grep -nosm 1 "using ParallelOption =[^;]*;$" ${1}.cpp)"
if [[ -z "$match" ]] ; then
    echo "$(basename "$0") : Parallelization option not found in ${1}.cpp. Declaration must be in a single line."
    exit 1
else
    lineno="$(cut -d ':' -f1 <<< "$match")"
    pattern="$(cut -d ':' -f 2- <<< "$match")"
    sed -i.bak "${lineno}s/${pattern}/using ParallelOption = par::ParallelOptions::${parallel};/" "${1}.cpp" && rm "${1}.cpp.bak"
    if [[ "$?" -ne 0 ]] ; then
	echo "$(basename "$0") : Failed to substitute parallel option in ${1}.cpp"
    fi
fi

match="$(grep -nosm 1 "using Model =[^;]*;$" ${1}.cpp)"
if [[ -z "$lineno" ]] ; then
    echo "$(basename "$0") : Model definition not found in ${1}.cpp. Declaration must be in a single line."
    exit 1
else
    lineno="$(cut -d ':' -f1 <<< "$match")"
    pattern="$(cut -d ':' -f 2- <<< "$match")"
    sed -i.bak "${lineno}s/${pattern}/using Model = ptof::Model::${2};/" "${1}.cpp" && rm "${1}.cpp.bak"
    if [[ "$?" -ne 0 ]] ; then
        echo "$(basename "$0") : Failed to substitute model definition in ${1}.cpp"
    fi
fi

exit 0
