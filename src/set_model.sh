#!/bin/bash
# Exit codes:
# 0 - Both models set successfuly
# 1 - Failed to set ParticleTrackingOF model, succeeded in setting ParticleTrackingOF_TwoPhaseStationary model
# 2 - Succeeded in setting ParticleTrackingOF model, failed to set ParticleTrackingOF_TwoPhaseStationary model
# 3 - Failed to set either model
# 4 - Bad parameters

if [[ "$#" -eq 1 ]] ; then
	err=0
	if ! grep -xsq ".*using namespace ptof::model_.*" ParticleTrackingOF.cpp ; then
		echo "set_model.sh : Model definition not found in ParticleTrackingOF.cpp. Using declaration must be in a single line."
		err=1
	fi
	sed -i'.bak' "s/.*using namespace ptof::model_.*/  using namespace ptof::model_$1;/" ParticleTrackingOF.cpp

	if ! grep -xsq ".*using namespace ptof::model_.*" ParticleTrackingOF_TwoPhaseNonStationary.cpp ; then
                echo "set_model.sh : Model definition not found in ParticleTrackingOF_TwoPhaseNonStationary.cpp. Using declaration must be in a single line."
		if [[ "$err" -eq 0 ]] ; then
			err=2
		else
			err=3
		fi
        fi
	sed -i'.bak' "s/.*using namespace ptof::model_.*/  using namespace ptof::model_$1;/" ParticleTrackingOF_TwoPhaseNonStationary.cpp
	rm *.bak

	exit $err
fi
echo "Bad parameters."
exit 4
