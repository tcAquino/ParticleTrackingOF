#!/bin/bash
if [[ "$#" -eq 1 ]]
then
	./set_model.sh $1
	make ParticleTrackingOF
	mv "../bin/ParticleTrackingOF" "../bin/ParticleTrackingOF_$1"
 ./set_model.sh $1
  make ParticleTrackingOF_TwoPhaseNonStationary
  mv "../bin/ParticleTrackingOF_TwoPhaseNonStationary" "../bin/ParticleTrackingOF_TwoPhaseNonStationary_$1"
	exit 0
fi
echo "Bad parameters."
exit 1
