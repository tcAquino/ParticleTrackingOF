GREEN='\033[1;32m'
RED='\033[1;31m'
NC='\033[0m'
err=0
./make_parallel.sh advection_diffusion_2d
((err+=$?))
./make_parallel.sh advection_diffusion_fpt_2d
((err+=$?))
./make_parallel.sh advection_diffusion_decay_catalytic_2d
((err+=$?))
./make_parallel.sh advection_2d
((err+=$?))
./make_parallel.sh periodic_cartesian_advection_diffusion_2d
((err+=$?))
./make_parallel.sh periodic_cartesian_advection_diffusion_fpt_2d
((err+=$?))
./make_parallel.sh periodic_cartesian_advection_diffusion_decay_catalytic_2d
((err+=$?))
./make_parallel.sh periodic_cartesian_advection_2d
((err+=$?))
./make_parallel.sh advection_diffusion_3d
((err+=$?))
./make_parallel.sh advection_diffusion_fpt_3d
((err+=$?))
./make_parallel.sh advection_diffusion_decay_catalytic_3d
((err+=$?))
./make_parallel.sh advection_3d
((err+=$?))
./make_parallel.sh periodic_cartesian_advection_diffusion_3d
((err+=$?))
./make_parallel.sh periodic_cartesian_advection_diffusion_fpt_3d
((err+=$?))
./make_parallel.sh periodic_cartesian_advection_diffusion_decay_catalytic_3d
((err+=$?))
./make_parallel.sh periodic_cartesian_advection_3d
((err+=$?))
./make_parallel.sh bcc_cartesian_advection_diffusion
((err+=$?))
./make_parallel.sh bcc_cartesian_advection_diffusion_fpt
((err+=$?))
./make_parallel.sh bcc_cartesian_advection_diffusion_decay_catalytic
((err+=$?))
./make_parallel.sh bcc_cartesian_advection
((err+=$?))
./make_parallel.sh bcc_symmetryplanes_advection_diffusion
((err+=$?))
./make_parallel.sh bcc_symmetryplanes_advection_diffusion_fpt
((err+=$?))
./make_parallel.sh bcc_symmetryplanes_advection_diffusion_decay_catalytic
((err+=$?))
./make_parallel.sh bcc_symmetryplanes_advection
((err+=$?))
if [[ err -ne 0 ]]
then
  echo -e "${RED}$err model compilations failed.${NC}"
  exit 1
else
  echo -e "${GREEN}All models compiled successfuly.${NC}"
  exit 0
fi
