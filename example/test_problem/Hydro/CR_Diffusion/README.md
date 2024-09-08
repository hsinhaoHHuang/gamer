# Compilation flags
- Must enable
   - [[MODEL=HYDRO | Installation: Simulation-Options#MODEL]]
   - [[EOS=EOS_GAMMA_CR | Installation: Simulation-Options#EOS]]
   - [[FLOAT8 | Installation: Simulation-Options#FLOAT8]]
   - [[COSMIC_RAY | Installation: Simulation-Options#COSMIC_RAY]]
   - [[CR_DIFFUSION | Installation: Simulation-Options#CR_DIFFUSION]]
- Must disable
   - [[COMOVING | Installation: Simulation-Options#COMOVING]]
   - [[PARTICLE | Installation: Simulation-Options#PARTICLE]]
   - [[GRAVITY | Installation: Simulation-Options#GRAVITY]]
- Available options
   - [[Miscellaneous Options | Installation: Simulation-Options#miscellaneous-options]]


# Default setup
1. `CR_Diffusion_Type`      4  (Blast wave)
2. `CR_Diffusion_Mag_Type`  4  (Suwa+ 2007)
3. 3D simulation

# Note
1. The analytical of `CR_Diffusion_Type={0,1,2,3}` need to set all fluid to be fixed but cosmic ray.


## `CR_Diffusion_Type`
1. Gaussian distribution ball (`CR_Diffusion_Type=0`)
   a. All of the fluid must be fixed, except for cosmic rays.

2. Step function ring         (`CR_Diffusion_Type=1`)
   a. All of the fluid must be fixed, except for cosmic rays.
   b. This test problem is a 2D simulation.
   c. The analytical solution assume the CR diffuse along azimuthal direction only.

3. Gaussian distribution ring (`CR_Diffusion_Type=2`)
   a. All of the fluid must be fixed, except for cosmic rays.
   b. This test problem is a 2D simulation.
   c. The analytical solution assumes the CRs diffuse along azimuthal direction only.

4. Gaussian distribution plane (`CR_Diffusion_Type=3`)
   a. All of the fluid must be fixed, except for cosmic rays.
   b. The analytical solution assumes the CRs diffuse along diagonal direction.

5. CR driven blast wave        (`CR_Diffusion_Type=4`)
   a. This test problem is a 3D simulation.


## `CR_Diffusion_Mag_Type`
1. Uniform
   a. Uniform magnetic field.
2. Circular
   a. Only works for 2D simulation
   b. x-y plane: `Bx=-y/r`, `By=x/r`, `r=sqrt(x^2+y^2)`
   c. Using `CR_Diffusion_MagX` as amplitude.
3. Random
   a. NOT divergence free!
   b. `-1. < B_i < 1.`
4. Radial
   a. NOT divergence free!
   b. Only works for 2D simulation
   c. x-y plane: `Bx=x/r`, `By=y/r`, `r=sqrt(x^2+y^2)`
   d. Using `CR_Diffusion_MagX` as amplitude.
5. Suwa+ 2007
   a. Follows Suwa+ 2007 setup
   b. Using `CR_Diffusion_MagX` as amplitude.


## Set simulation dimensions
The dimensions are controlled by `CR_Diffusion_G{X,Y,Z}`.

For example, the 1D simulation running on y-axis should set the followings:
  ```
  CR_diffusion_GX 0
  CR_diffusion_GY 1
  CR_diffusion_GZ 0
  ```

For example, the 2D simulation running on x-z plane should set the followings:
  ```
  CR_diffusion_GX 1
  CR_diffusion_GY 0
  CR_diffusion_GZ 1
  ```


## How to fix fluid
Not available right now.


## `yt_L1error`
1. The file structure should be like:
   <pre>
   ./ --+-- res_0032 --+-- Data_000000
        |              |
        |              |-- ...
        |
        |-- res_0064 --+-- Data_000000
        |              |
        |-- ...        |-- ...
   <\pre>
2. Uncomment or comment your simulation type and set the parameters.
