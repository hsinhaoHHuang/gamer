# Compilation flags
- Must enable
   - [[MODEL=ELBDM | Installation: Simulation-Options#MODEL]]
   - [[GRAVITY | Installation: Simulation-Options#GRAVITY]]
   - [[COMOVING | Installation: Simulation-Options#COMOVING]]
- Must disable
   - [[PARTICLE | Installation: Simulation-Options#PARTICLE]]
   - [[UNSPLIT_GRAVITY | Installation: Simulation-Options#UNSPLIT_GRAVITY]]
- Available options
   - [[Miscellaneous Options | Installation: Simulation-Options#miscellaneous-options]]


# Default setup
1. 32^3 uniform resolution
   --> No grid refinement is adopted in this test

2. [[DT__MAX_DELTA_A | Runtime-Parameters:-Timestep#DT__MAX_DELTA_A]] is set to an extremely large number so that the time-step will
   be controlled by the fluid solver

3. The default setup in `Input__Parameter` and `Input__TestProb` gives an UNSTABLE
   solution, in which the solution grows exponentially. To test a STABLE solution,
   one can try the following parameters:

   * `Input__Parameter`:
      - BOX_SIZE  = 0.5
      - A_INIT    = 1.0e-5
      - OUTPUT_DT = 2.0e-6

   * `Input__TestProb`:
      - Jeans_RealAmp0 = 4.0e-8

> [!CAUTION]
> This setup requires further verification


# Note
1. Analytical solution reference: [Woo, T. & Chiueh, T. 2009, ApJ, 697, 850](https://doi.org/10.1088/0004-637X/697/1/850)

> [!NOTE]
> The imaginary part `I` grows much faster than the real part `R`.
  Consequently, the accuracy of linear prediction will be deteriorated
  when the assumption (`2R >> I^2`) starts to break down, especially in the
  higher-resolution tests, and hence the 2nd-order accuracy may no longer
  hold.

> [!CAUTION]
> 2nd-order accuracy has NOT been verified for the UNSTABLE solution