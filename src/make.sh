#!/bin/bash
if [[ "$#" -eq 1 ]]
then
	./set_model.sh $1
	make ParticleTrackingOF
	mv "../bin/ParticleTrackingOF" "../bin/ParticleTrackingOF_$1"
	exit 0
fi
echo "Bad parameters."
exit 1
