#!/bin/bash
if [[ "$#" -eq 1 ]]
then
	sed -i'.bak' "s/.*using namespace ptof::model_.*/  using namespace ptof::model_$1;/" ParticleTrackingOF.cpp
	rm *.bak
	exit 0
fi
echo "Bad parameters."
exit 1
