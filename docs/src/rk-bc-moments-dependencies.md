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
2. In moment-kinetic cases, update the electron moments that are calculated
   directly from ion moments, `scratch.electron_density` and
   `scratch.electron_n_upar`.
3. Calculate `moments.ion.upar`, `moments.ion.vth`, `moments.electron.upar`,
   `moments.electron.vth`, moments.neutral.uz`, and `moments.neutral.vth` from
   evolved moments.
   - `moments.electron.upar` depends on
     `scratch.electron_density/electron_n_upar`, and so can only be updated
     here in moment-kinetic cases.
4. Set negative values in the shape functions to zero.
5. Apply sheath entrance boundary condition for electrons.
   - Will eventually have to calculate `phi` and `this_scratch.electron_n_upar`
     from electron parallel force balance, consistently with parallel currents
     through sheaths and `phi` at sheath entrances. Will then require
     electron-ion and electron-neutral friction, so likely to depend on
     `moments.ion.upar` and `moments.neutral.uz`, but (if including $E\times B$
     drifts) also determines the ion flux to the wall, and therefore the
     neutral flux from the wall. This is therefore likely to require an
     iteration to determine a consistent set of fluxes and phi.
   - as a first step, before the iteration in the previous point is
     implemented, apply the boundary condition on `scratch.neutral_n_uz` and
     `moments.neutral.uz` before the bc, but the boundary condition on the
     neutral shape function (which would need to be done after the iteration)
     after this bc.
   - requires negative values in electron shape function to be set to zero
     first.
   - kinetic electrons require fully moment-kinetic version, so this bc does
     not depend on the ion or neutral shape functions (or their boundary
     conditions).
6. Apply moment constraint corrections to electron shape function.
   - must be done after applying boundary conditions to electron shape function
     (the boundary conditions take care of enforcing constraints at the
     boundary points, and must not modify the shape functions on the boundary
     points before applying boundary conditions).
7. (Unless full-f) calculate `phi`, `Ez`, vEr`, `vEz`.
   - depends on `scratch.electron_density`.
   - eventually `phi` calculated as part of iteration with electron shape
     function boundary conditions, when using kinetic electrons, instead of
     here.
8. Derived moments for electrons, e.g. `qpar` and `ppar`/`pperp`. Derivatives of
   moments.
   - depends on fully updated electron shape functions and input moments,
     including any moment constraints and boundary conditions.
9. Apply the sheath-entrance boundary condition for ions.
   - sets ion shape function at sheath entrance when ions are treated
     kinetically, or `this_scratch.n_upar` and `moments.ion.upar` at the sheath
     entrance when ions are treated as fluid.
   - depends on `moments.ion.upar` and `moments.ion.vth` for kinetic version,
     but only RK updated evolving variables for fluid version.
   - depends on `phi` when including $E\times B$ drifts. For full-f version
     (with Boltzmann electron response) this makes boundary condition poorly
     defined, because phi is needed to apply bc on $f_i$, but phi depends on
     density and density is (possibly) modified by applying the boundary
     condition. We fudge this issue by using the old phi to apply the boundary
     condition, so that at least density and phi are consistent with $f_i$ at
     the end of each RK stage, even though the boundary condition was not
     applied fully consistently.
   - requires negative values in ion shape function to be set to zero first.
10. Apply moment constraint corrections to ion shape function.
   - must be done after applying boundary conditions to ion shape function (the
     boundary conditions take care of enforcing constraints at the boundary
     points, and must not modify the shape functions on the boundary points
     before applying boundary conditions).
11. Derived moments for ions, e.g. `qpar` and `ppar`/`pperp`. Derivatives of
   moments. In `evolve_n_upar=false` case, this includes calculation of the ion
   parallel particle flux.
   - depends on fully updated ion shape functions and input moments, including
     any moment constraints and boundary conditions.
12. Apply the wall boundary condition for neutrals. Sets boundary values of
   `this_scratch.n_uz_neutral` and `moments.neutral.uz`.
   - depends on ion particle flux to wall (in non-moment-kinetic case, this
     moment must have been calculated after applying the ion bc), and
     `moments.neutral.vth` (in the `evolve_p=true` case).
   - requires negative values in neutral shape function to be set to zero
     first.
13. Apply moment constraint corrections to neutral shape function.
   - must be done after applying boundary conditions to neutral shape function
     (the boundary conditions take care of enforcing constraints at the
     boundary points, and must not modify the shape functions on the boundary
     points before applying boundary conditions).
14. Derived moments for neutrals, e.g. `qz` and `pz`/`pzeta`/`pr`. Derivatives
   of moments.
   - depends on fully updated neutral shape function and input moments,
     including any moment constraints and boundary conditions.
15. When using kinetic electrons, fill the `moments` struct with the RK stage
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
