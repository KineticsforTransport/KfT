# Dependencies between RK updates, boundary conditions and moment calculations

It is important that updates are done in the correct order to ensure
consistency, otherwise the accuracy of the Runge-Kutta (RK) timestep is
degraded to first order. This page documents the dependencies and the intended
order of operations so that the source code can be checked, considering the
fully moment-kinetic model (the 'full-f' model is simpler as there are fewer
dependencies on moments).

1. RK updates of evolving variables (at grid points that are not set by
   boundary conditions).
   - no dependencies.
   - all other steps require this to be done first.
2. Calculate `moments.ion.upar`, `moments.ion.vth`, `moments.electron.upar`,
   `moments.electron.vth`, moments.neutral.uz`, and `moments.neutral.vth` from
   evolved moments.
3. Set negative values in the shape functions to zero.
4. Apply the sheath-entrance boundary condition for ions.
   - sets ion shape function at sheath entrance when ions are treated
     kinetically, or `this_scratch.n_upar` and `moments.ion.upar` at the sheath
     entrance when ions are treated as fluid.
   - depends on `moments.ion.upar` and `moments.ion.vth` for kinetic version,
     but only RK updated evolving variables for fluid version.
   - requires negative values in ion shape function to be set to zero first.
5. Apply the wall boundary condition for neutrals. Sets
   `this_scratch.n_uz_neutral` and `moments.neutral.uz`.
   - depends on ion particle flux to wall, and `moments.neutral.vth`.
   - requires negative values in neutral shape function to be set to zero
     first.
6. Apply sheath entrance boundary condition for electrons.
   - Will eventually have to calculate `phi` and `this_scratch.electron_n_upar`
     from electron parallel force balance, consistently with parallel currents
     through sheaths and `phi` at sheath entrances. Will then require
     electron-ion and electron-neutral friction, so likely to depend on
     `moments.ion.upar` and `moments.neutral.uz`.
   - requires negative values in electron shape function to be set to zero
     first.
7. Apply moment constraint corrections to shape functions.
   - must be done after applying boundary conditions to shape functions (the
     boundary conditions take care of enforcing constraints at the boundary
     points, and must not modify the shape functions on the boundary points
     before applying boundary conditions).
8. Other derived moments, e.g. `qpar` and `ppar`/`pperp`. Derivatives of
   moments.
   - depends on fully updated shape functions and input moments, including any
     moment constraints and boundary conditions.
9. When using kinetic electrons, fill the `moments` struct with the RK stage
   values of the evolving variables. This is needed as to be compatible with
   both pseudo-transient (which has a special `scratch_electron` struct just
   for electrons) and time-evolving approaches (which use the `scratch` struct
   from the main RK solve), the electron implicit solves reads several
   variables from the `moments` struct, which would not otherwise be updated
   until the end of the RK timestep (rather than at each stage). The variables
   that need to be (or might need to be) updated are `moments.ion.dens`,
   `moments.ion.n_upar`, `moments.ion.p`, `moments.electron.dens`,
   `moments.electron.upar`, `moments.electron.p`, `moments.neutral.dens`,
   `moments.neutral.n_uz`, `moments.neutral.p`.
