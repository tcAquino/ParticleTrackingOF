#!/bin/bash
if [[ "$#" -eq 1 ]]
then
	sed -i'.bak' "s/.*using namespace ptof::model_.*/  using namespace ptof::model_$1_parallel;/" ParticleTrackingOF_Parallel.cpp
  sed -i'.bak' "s/.*using namespace ptof::model_.*/  using namespace ptof::model_$1_parallel;/" ParticleTrackingOF_TwoPhaseNonStationary_Parallel.cpp
	rm *.bak
	exit 0
fi
echo "Bad parameters."
exit 1
