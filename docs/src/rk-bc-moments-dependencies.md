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
3. Apply the sheath-entrance boundary condition for ions.
   - sets ion shape function at sheath entrance when ions are treated
     kinetically, or `this_scratch.n_upar` and `moments.ion.upar` at the sheath
     entrance when ions are treated as fluid.
   - depends on `moments.ion.upar` and `moments.ion.vth` for kinetic version,
     but only RK updated evolving variables for fluid version.
4. Apply the wall boundary condition for neutrals. Sets
   `this_scratch.n_uz_neutral` and `moments.neutral.uz`.
   - depends on ion particle flux to wall, and `moments.neutral.vth`.
5. Apply sheath entrance boundary condition for electrons.
   - Will eventually have to calculate `phi` and `this_scratch.electron_n_upar`
     from electron parallel force balance, consistently with parallel currents
     through sheaths and `phi` at sheath entrances. Will then require
     electron-ion and electron-neutral friction, so likely to depend on
     `moments.ion.upar` and `moments.neutral.uz`.
6. Apply moment constraint corrections to shape functions.
   - must be done after applying boundary conditions to shape functions (the
     boundary conditions take care of enforcing constraints at the boundary
     points, and must not modify the shape functions on the boundary points
     before applying boundary conditions).
7. Other derived moments, e.g. `qpar` and `ppar`/`pperp`. Derivatives of
   moments.
   - depends on fully updated shape functions and input moments, including any
     moment constraints and boundary conditions.
