# This script should run in the same directory as configure.py

PYTHON=python

${PYTHON} configure.py --machine=eureka_intel --model=HYDRO --flu_scheme=MHM_RP --flux=HLLD --mhd=true \
                       --microphysics=true --cosmic_ray=true --eos=COSMIC_RAY --cr_diffusion=true
