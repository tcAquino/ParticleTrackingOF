#!/bin/bash
if [[ "$#" -eq 1 ]]
then
	./set_model_parallel.sh $1
	make ParticleTrackingOF_Parallel
	mv "../bin/ParticleTrackingOF_Parallel" "../bin/ParticleTrackingOF_$1_parallel"
 ./set_model_parallel.sh $1
  make ParticleTrackingOF_TwoPhaseNonStationary_Parallel
  mv "../bin/ParticleTrackingOF_TwoPhaseNonStationary" "../bin/ParticleTrackingOF_TwoPhaseNonStationary_$1_parallel"
	exit 0
fi
echo "Bad parameters."
exit 1
