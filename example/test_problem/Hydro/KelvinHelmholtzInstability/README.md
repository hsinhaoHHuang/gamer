# `configure.py` options
- Must enable
   - [[--model | Installation:-Option-List#--model]]=`HYDRO`
- Must disable
   - [[--comoving | Installation:-Option-List#--comoving]]
   - [[--particle | Installation:-Option-List#--particle]]
   - [[--gravity | Installation:-Option-List#--gravity]]
- Available options
   - [[Miscellaneous Options | Installation:-Option-List#miscellaneous-options]]


# Default setup
1. Refinement criteria: vorticity and regions near the shear plane
2. Maximum refinement level ([[MAX_LEVEL | Runtime-Parameters:-Refinement#MAX_LEVEL]])=`3`


# Note
1. Adopt the periodic boundary condition
2. Shear velocity is defined on the x-y plane
3. Random values are added to velocity in all three directions
4. yt script `plot_density.py` for visualization
5. Must disable `OPT__INIT_GRID_WITH_OMP`
   - Otherwise all threads would share the same random seed
