module implicit_ion_parallel_dynamics_functions

export ion_parallel_dynamics_backward_euler!

using ..alpha_advection: update_speed_alpha_no_sr!
using ..boundary_conditions
using ..calculus: reconcile_element_boundaries_MPI_z_pdf_vpavperpz!
using ..charge_exchange: charge_exchange_collisions_mk_single_species_r!
using ..communication
using ..continuity: continuity_equation_no_sr!
using ..energy_equation: energy_equation_no_sr!
using ..external_sources: total_external_Maxwellian_sources_no_sr!
using ..force_balance: force_balance_no_sr!
using ..ion_jacobian_terms
using ..ionization: ionization_collisions_mk_ion_no_sr!
using ..jacobian_matrices: get_sparse_joined_array
using ..krook_collisions: krook_collisions_no_sr!
using ..looping
using ..moment_kinetics_structs
using ..nonlinear_solvers
using ..numerical_dissipation: vpa_dissipation_no_sr!, vperp_dissipation_no_sr!
using ..r_advection: update_speed_r_no_sr!
using ..source_terms: source_terms_no_sr!
using ..timer_utils
using ..type_definitions
using ..velocity_moments: update_ppar_species_no_sr!, calculate_ion_qpar_from_pdf_no_r!,
                          calculate_implicit_ion_moment_derivatives!,
                          update_derived_ion_moment_time_derivatives_no_sr!
using ..vpa_advection: vpa_advection_no_sr!
using ..vperp_advection: vperp_advection_no_sr!
using ..z_advection: z_advection_no_sr!

using ConcreteStructs
using LinearAlgebra

"""
    ion_parallel_dynamics_backward_euler!(
        fvec_out, fvec_in, pdf, moments, fields, composition, collisions,
        geometry, external_source_settings, manufactured_source_list,
        fp_arrays, t_params, scratch_dummy, r::coordinate, z::coordinate,
        vperp::coordinate, vpa::coordinate, vzeta::coordinate,
        vr::coordinate, vz::coordinate, gyrophase::coordinate,
        spectral_objects, advection_structs, num_diss_params,
        advance_implicit, nl_solver_params, dt, is, ir)

Take a single backward euler timestep for the ion shape function \$g_i\$.
"""
@timeit global_timer ion_parallel_dynamics_backward_euler!(
                         fvec_out, fvec_in, pdf, moments, fields, composition, collisions,
                         geometry, external_source_settings, manufactured_source_list,
                         fp_arrays, t_params, scratch_dummy, r::coordinate, z::coordinate,
                         vperp::coordinate, vpa::coordinate, vzeta::coordinate,
                         vr::coordinate, vz::coordinate, gyrophase::coordinate,
                         spectral_objects, advection_structs, num_diss_params,
                         advance_implicit, nl_solver_params, dt, is, ir) = begin

    # create several 0D dummy arrays for use in taking derivatives
    buffer_1 = @view scratch_dummy.buffer_rs_1[ir,is]
    buffer_2 = @view scratch_dummy.buffer_rs_2[ir,is]
    buffer_3 = @view scratch_dummy.buffer_rs_3[ir,is]
    buffer_4 = @view scratch_dummy.buffer_rs_4[ir,is]

    f_old = @view fvec_in.pdf[:,:,:,ir,is]
    dens = @view fvec_in.density[:,ir,is]
    upar = @view fvec_in.upar[:,ir,is]
    p = @view fvec_in.p[:,ir,is]
    dens_neutral = @view fvec_in.density_neutral[:,ir,is]
    uz_neutral = @view fvec_in.uz_neutral[:,ir,is]
    p_neutral = @view fvec_in.p_neutral[:,ir,is]
    f_new = @view fvec_out.pdf[:,:,:,ir,is]

    left_preconditioner, right_preconditioner, recalculate_preconditioner! =
        @views get_ion_parallel_dynamics_preconditioner(
            nl_solver_params, f_new, buffer_1, buffer_2, buffer_3, buffer_4, dens, upar,
            p, dens_neutral, uz_neutral, p_neutral, moments, fields, composition,
            collisions, geometry, z, vperp, vpa, spectral_objects.z_spectral,
            spectral_objects.vperp_spectral, spectral_objects.vpa_spectral,
            advection_structs.r_advect[:,:,:,ir,is],
            advection_structs.alpha_advect[:,:,:,ir,is],
            advection_structs.z_advect[:,:,:,ir,is],
            advection_structs.vperp_advect[:,:,:,ir,is],
            advection_structs.vpa_advect[:,:,:,ir,is], scratch_dummy,
            external_source_settings, num_diss_params, t_params, is, ir)

    # Does preconditioner need to be recalculated?
    if t_params.dt[] > 1.5 * nl_solver_params.precon_dt[] ||
            t_params.dt[] < 2.0/3.0 * nl_solver_params.precon_dt[]

        # dt has changed significantly, so update the preconditioner
        nl_solver_params.solves_since_precon_update[] = nl_solver_params.preconditioner_update_interval
    end
    if nl_solver_params.solves_since_precon_update[] ≥ nl_solver_params.preconditioner_update_interval
        recalculate_preconditioner!()
    end

    # Do a backward-Euler update of the ion shape function.
    residual_func! = ion_parallel_dynamics_residual!(
                         fvec_in, pdf, moments, fields, composition, collisions,
                         geometry, external_source_settings, manufactured_source_list,
                         fp_arrays, t_params, scratch_dummy, r, z, vperp, vpa, vzeta, vr,
                         vz, gyrophase, spectral_objects, advection_structs,
                         num_diss_params, advance_implicit,
                         nl_solver_params.linear_iterations, dt, is, ir)

    residual = scratch_dummy.implicit_buffer_vpavperpz_1
    delta_x = scratch_dummy.implicit_buffer_vpavperpz_2
    rhs_delta = scratch_dummy.implicit_buffer_vpavperpz_3
    v = scratch_dummy.implicit_buffer_vpavperpz_4
    w = scratch_dummy.implicit_buffer_vpavperpz_5

    newton_success = newton_solve!(f_new, residual_func!, residual, delta_x, rhs_delta, v,
                                   w, nl_solver_params;
                                   left_preconditioner=left_preconditioner,
                                   right_preconditioner=right_preconditioner,
                                   coords=(z=z, vperp=vperp, vpa=vpa))

    if !newton_success
        # Newton solve failed, so update the preconditioner
        nl_solver_params.solves_since_precon_update[] = nl_solver_params.preconditioner_update_interval
    end

    return newton_success
end

"""
ion_parallel_dynamics_euler_time_advance_no_sr! advances the vector
equation dfvec/dt = G[f] that includes the z and velocity space terms
in the ion kinetic equation using the forward Euler method:
fvec_out = fvec_in + dt*fvec_in,
with fvec_in an input and fvec_out the output.

Note `dt` is passed separately from `t_params` because sometimes (in the IMEX Runge-Kutta
implementation), a call needs to be made with `dt` scaled by some coefficient.
"""
@timeit global_timer ion_parallel_dynamics_euler_time_advance_no_sr!(
                         fvec_out, fvec_in, pdf, debug_scratch_electron, fields, moments,
                         advect_objects, vz, vr, vzeta, vpa, vperp, gyrophase, z, r,
                         t_params, dt, spectral_objects, composition, collisions,
                         geometry, scratch_dummy, manufactured_source_list,
                         external_source_settings, num_diss_params, advance, fp_arrays,
                         solve_counter, is, ir) = begin

    # Convenience wrapper for calls to write debug information within this function.
    # Not sure the debug-IO functionality could be used in this function if r.n>1, so
    # comment out for now. If running in 1D (r.n==1), can uncomment
    # write_debug_data_to_binary() below and debug-IO stuff should work...
    function write_debug_IO(label)
        if t_params.debug_io === nothing
            # Allow compiler to optimise away this function if debug IO is not being used.
            return nothing
        end
        #write_debug_data_to_binary(fvec_out, debug_scratch_electron, moments, fields,
        #                           composition, t_params, r, z, vperp, vpa, vzeta, vr, vz,
        #                           label, solve_counter)
        return nothing
    end
    write_debug_IO("begin ion_parallel_dynamics_euler_time_advance_no_sr!")

    # define some abbreviated variables for tidiness
    n_ion_species = composition.n_ion_species
    n_neutral_species = composition.n_neutral_species
    t = t_params.t[]

    vpa_spectral, vperp_spectral, z_spectral = spectral_objects.vpa_spectral, spectral_objects.vperp_spectral, spectral_objects.z_spectral
    vz_spectral = spectral_objects.vz_spectral
    vpa_advect = @view advect_objects.vpa_advect[:,:,:,ir,is]
    vperp_advect = @view advect_objects.vperp_advect[:,:,:,ir,is]
    r_advect = @view advect_objects.r_advect[:,:,:,ir,is]
    alpha_advect = @view advect_objects.alpha_advect[:,:,:,ir,is]
    z_advect = @view advect_objects.z_advect[:,:,:,ir,is]

    f_in = @view fvec_in.pdf[:,:,:,ir,is]
    dens = @view fvec_in.density[:,ir,is]
    upar = @view fvec_in.upar[:,ir,is]
    p = @view fvec_in.p[:,ir,is]
    vth = @view moments.ion.vth[:,ir,is]
    ppar = @view moments.ion.ppar[:,ir,is]
    qpar = @view moments.ion.qpar[:,ir,is]
    evolve_density = moments.evolve_density
    evolve_upar = moments.evolve_upar
    evolve_p = moments.evolve_p

    # Update the ppar and qpar moments that depend on the shape function, and the moment
    # time derivatives that will be used in the shape function equation.
    update_ppar_species_no_sr!(ppar, dens, upar, vth, p, f_in, vpa, vperp, z,
                               evolve_density, evolve_upar, evolve_p)
    calculate_ion_qpar_from_pdf_no_r!(qpar, dens, upar, vth, f_in, vpa, vperp, z,
                                      evolve_density, evolve_upar, evolve_p)
    calculate_implicit_ion_moment_derivatives!(moments, scratch_dummy, z, vperp,
                                               z_spectral, is, ir)
    continuity_equation_no_sr!(fvec_in, fields, moments, composition, geometry, dt,
                               collisions.reactions.ionization_frequency,
                               external_source_settings.ion, num_diss_params, is, ir)
    force_balance_no_sr!(fvec_in, moments, fields, collisions, dt, composition, geometry,
                         external_source_settings.ion, num_diss_params, z, is, ir)
    energy_equation_no_sr!(fvec_in, moments, fields, collisions, dt, composition,
                           geometry, external_source_settings.ion, num_diss_params, is,
                           ir)
    update_derived_ion_moment_time_derivatives_no_sr!(fvec_in, moments, is, ir)

    # z_advection! advances 1D advection equation in z
    # apply z-advection operation to ion species
    if advance.z_advection
        z_advection_no_sr!(fvec_out.pdf, fvec_in, moments, fields, z_advect, z, vpa,
                           vperp, dt, z_spectral, composition, geometry, scratch_dummy,
                           is, ir)
        write_debug_IO("ion_parallel_dynamics z_advection_no_sr!")
    end

    if advance.vpa_advection || advance.vperp_advection
        # Update alpha and r advection speeds, that are needed for vpa_advection.
        update_speed_alpha_no_sr!(alpha_advect, evolve_upar, evolve_p, fields, vpa, vperp,
                                  z, geometry, is, ir)
        update_speed_r_no_sr!(r_advect, fields, evolve_density, evolve_upar, evolve_p,
                              vpa, vperp, r, z, geometry, is, ir)
    end

    # vpa_advection! advances the 1D advection equation in vpa. Note this must be called
    # after z_advection!() so that the z-advection speed has been updated, and after
    # updating alpha and r advection speeds.
    if advance.vpa_advection
        vpa_advection_no_sr!(fvec_out.pdf, fvec_in, fields, moments, vpa_advect, r_advect,
                             alpha_advect, z_advect, vpa, vperp, z, dt, t,
                             vpa_spectral, composition, collisions,
                             external_source_settings.ion, geometry, is, ir)
        write_debug_IO("ion_parallel_dynamics vpa_advection_no_sr!")
    end

    # vperp_advection requires information about z, alpha, and r advection
    # so call vperp_advection! only after z advection and alpha, r speed routines
    if advance.vperp_advection
        vperp_advection_no_sr!(fvec_out.pdf, fvec_in, vperp_advect, r, z, vperp, vpa, dt,
                               vperp_spectral, composition, z_advect, alpha_advect,
                               r_advect, geometry, moments, fields, is, ir)
        write_debug_IO("ion_parallel_dynamics vperp_advection_no_sr!")
    end

    if advance.source_terms
        source_terms_no_sr!(fvec_out.pdf, fvec_in, moments, r_advect, alpha_advect,
                            z_advect, vpa, vperp, z, dt, composition, collisions,
                            external_source_settings.ion, is, ir)
        write_debug_IO("ion_parallel_dynamics source_terms_no_sr!")
    end

    if advance.ion_cx_collisions || advance.ion_ionization_collisions
        error("2V/3V ion/neutral collsions not currently supported in "
              * "ion_parallel_dynamics_euler_time_advance_no_sr!()")
        # gyroaverage neutral dfn and place it in the ion.buffer array for use in the collisions step
        vzvrvzeta_to_vpavperp_no_sr!(pdf.ion.buffer, fvec_in.pdf_neutral, vz, vr, vzeta,
                                     vpa, vperp, gyrophase, z, geometry, composition, is,
                                     ir)
        write_debug_IO("ion_parallel_dynamics vzvrvzeta_to_vpavperp_no_sr!")
    end

    # account for charge exchange collisions between ions and neutrals
    if advance.ion_cx_collisions_1V
        if !evolve_density
            error("charge_exchange_collisions_mk_single_species_r!() is not the right "
                  * "method for 'full-f' mode with evolve_density=false. Don't think we "
                  * "need full-f support here, so unsupported at least for now.")
        end
        @begin_anyzv_z_region()
        @views charge_exchange_collisions_mk_single_species_r!(
                   fvec_out.pdf[:,1,:], fvec_in.pdf[:,1,:], fvec_in.pdf_neutral[:,1,1,:],
                   fvec_in.density_neutral, upar, fvec_in.uz_neutral, vth,
                   moments.neutral.vth[:,ir,is], moments, vpa, vz,
                   collisions.reactions.charge_exchange_frequency, vz_spectral, dt;
                   neutrals=false)
        write_debug_IO("ion_parallel_dynamics charge_exchange_collisions_mk_single_species_r!")
    elseif advance.ion_cx_collisions
        error("2V/3V ion/neutral collsions not currently supported in "
              * "ion_parallel_dynamics_euler_time_advance_no_sr!()")
        ion_charge_exchange_collisions_3V_no_sr!(fvec_out.pdf, pdf.ion.buffer, fvec_in,
                                                 composition, vz, vr, vzeta, vpa, vperp,
                                                 z,
                                                 collisions.reactions.charge_exchange_frequency,
                                                 dt, is, ir)
        write_debug_IO("ion_parallel_dynamics_ion_charge_exchange_collisions_3V!")
    end
    # account for ionization collisions between ions and neutrals
    if advance.ion_ionization_collisions_1V
        @begin_anyzv_z_region()
        @views ionization_collisions_mk_ion_no_sr!(
                   fvec_out.pdf[:,1,:], fvec_in.pdf[:,1,:], fvec_in.pdf_neutral[:,1,1,:],
                   fvec_in.density_neutral, upar, fvec_in.uz_neutral,
                   collisions.reactions.ionization_frequency, moments, vpa, vz,
                   vz_spectral, dt, is, ir)
        write_debug_IO("ion_parallel_dynamics ion_ionization_collisions_1V_no_sr!")
    elseif advance.ion_ionization_collisions
        error("2V/3V ion/neutral collsions not currently supported in "
              * "ion_parallel_dynamics_euler_time_advance_no_sr!()")
        ion_ionization_collisions_3V_no_sr!(fvec_out.pdf, pdf.ion.buffer, fvec_in,
                                            composition, vz, vr, vzeta, vpa, vperp, z,
                                            collisions, dt, is, ir)
        write_debug_IO("ion_parallel_dynamics_ion_ionization_collisions_3V!")
    end

    # Add Krook collision operator for ions
    if advance.krook_collisions_ii
        krook_collisions_no_sr!(fvec_out.pdf, fvec_in, moments, composition, collisions,
                                vperp, vpa, dt, is, ir)
        write_debug_IO("ion_parallel_dynamics krook_collisions_no_sr!")
    end

    if advance.external_source
        total_external_Maxwellian_sources_no_sr!(fvec_out.pdf, fvec_in, moments,
                                                 external_source_settings.ion, vperp, vpa,
                                                 dt, is, ir)
        write_debug_IO("ion_parallel_dynamics_total_external_ion_sources!")
    end

    # add numerical dissipation
    if advance.ion_numerical_dissipation_vpa
        vpa_dissipation_no_sr!(fvec_out.pdf, fvec_in.pdf, vpa, vpa_spectral, dt,
                               num_diss_params.ion.vpa_dissipation_coefficient)
        write_debug_IO("ion_parallel_dynamics vpa_dissipation_no_sr!")
    end
    if advance.ion_numerical_dissipation_vperp
        vperp_dissipation_no_sr!(fvec_out.pdf, fvec_in.pdf, vperp, vperp_spectral, dt,
                                 num_diss_params.ion.vperp_dissipation_coefficient)
        write_debug_IO("ion_parallel_dynamics vperp_dissipation_no_sr!")
    end
    if advance.ion_numerical_dissipation_z
        error("z_dissipation not currently supported in "
              * "ion_parallel_dynamics_euler_time_advance_no_sr!()")
        z_dissipation_no_sr!(fvec_out.pdf, fvec_in.pdf, z, z_spectral, dt,
                             num_diss_params.ion.z_dissipation_coefficient,
                             num_diss_params.ion.z_dissipation_degree, scratch_dummy)
        write_debug_IO("ion_parallel_dynamics z_dissipation_no_sr!")
    end

    return nothing
end

@kwdef @concrete struct ion_parallel_dynamics_recalculate_lu!
    nl_solver_params
    f_new
    dens
    upar
    p
    dens_neutral
    uz_neutral
    p_neutral
    moments
    fields
    composition
    collisions
    geometry
    z
    vperp
    vpa
    z_spectral
    vperp_spectral
    vpa_spectral
    r_advect
    alpha_advect
    z_advect
    vperp_advect
    vpa_advect
    scratch_dummy
    external_source_settings
    num_diss_params
    t_params
    is::mk_int
    ir::mk_int
end

function (recalc::ion_parallel_dynamics_recalculate_lu!)()
    nl_solver_params = recalc.nl_solver_params
    f_new = recalc.f_new
    dens = recalc.dens
    upar = recalc.upar
    p = recalc.p
    dens_neutral = recalc.dens_neutral
    uz_neutral = recalc.uz_neutral
    p_neutral = recalc.p_neutral
    moments = recalc.moments
    fields = recalc.fields
    composition = recalc.composition
    collisions = recalc.collisions
    geometry = recalc.geometry
    z = recalc.z
    vperp = recalc.vperp
    vpa = recalc.vpa
    z_spectral = recalc.z_spectral
    vperp_spectral = recalc.vperp_spectral
    vpa_spectral = recalc.vpa_spectral
    r_advect = recalc.r_advect
    alpha_advect = recalc.alpha_advect
    z_advect = recalc.z_advect
    vperp_advect = recalc.vperp_advect
    vpa_advect = recalc.vpa_advect
    scratch_dummy = recalc.scratch_dummy
    external_source_settings = recalc.external_source_settings
    num_diss_params = recalc.num_diss_params
    t_params = recalc.t_params
    is = recalc.is
    ir = recalc.ir

global_rank[] == 0 && println("recalculating ion precon")

    dt = t_params.dt[]
    nl_solver_params.solves_since_precon_update[] = 0
    nl_solver_params.precon_dt[] = dt

    orig_lu, precon, input_buffer, output_buffer =
        nl_solver_params.preconditioners[ir]

    fill_ion_kinetic_equation_Jacobian_evolve_nup!(
        precon, f_new, dens, upar, p, dens_neutral, uz_neutral, p_neutral, moments,
        fields, composition, collisions, external_source_settings, geometry,
        num_diss_params, scratch_dummy, r_advect, alpha_advect, z_advect, vperp_advect,
        vpa_advect, z, vperp, vpa, z_spectral, vperp_spectral, vpa_spectral, dt, is, ir)

    @begin_anyzv_region()
    if anyzv_subblock_rank[] == 0
        if size(orig_lu) == (1, 1)
            # Have not properly created the LU decomposition before, so
            # cannot reuse it.
            @timeit_debug global_timer "lu" nl_solver_params.preconditioners[ir] =
                (lu(get_sparse_joined_array(precon)), precon, input_buffer, output_buffer)
        else
            # LU decomposition was previously created. The Jacobian always
            # has the same sparsity pattern, so by using `lu!()` we can
            # reuse some setup.
            try
                @timeit_debug global_timer "lu!" lu!(orig_lu, get_sparse_joined_array(precon); check=false)
            catch e
                if !isa(e, ArgumentError)
                    rethrow(e)
                end
                println("Sparsity pattern of matrix changed, rebuilding "
                        * "LU from scratch")
                @timeit_debug global_timer "lu" orig_lu = lu(get_sparse_joined_array(precon))
            end
            nl_solver_params.preconditioners[ir] = (orig_lu, precon, input_buffer,
                                                    output_buffer)
        end
    else
        nl_solver_params.preconditioners[ir] = (orig_lu, precon, input_buffer,
                                                output_buffer)
    end

    return nothing
end

@kwdef @concrete struct ion_parallel_dynamics_lu_precon!
    nl_solver_params
    scratch_dummy
    buffer_1
    buffer_2
    buffer_3
    buffer_4
    z
    is::mk_int
    ir::mk_int
end

@timeit_debug global_timer (lu_pre::ion_parallel_dynamics_lu_precon!)(precon_f) = begin
    nl_solver_params = lu_pre.nl_solver_params
    scratch_dummy = lu_pre.scratch_dummy
    buffer_1 = lu_pre.buffer_1
    buffer_2 = lu_pre.buffer_2
    buffer_3 = lu_pre.buffer_3
    buffer_4 = lu_pre.buffer_4
    z = lu_pre.z
    is = lu_pre.is
    ir = lu_pre.ir

    precon_lu, _, this_input_buffer, this_output_buffer =
        nl_solver_params.preconditioners[ir,is]

    @begin_anyzv_region()
    counter = 1
    @loop_z_vperp_vpa iz ivperp ivpa begin
        this_input_buffer[counter] = precon_f[ivpa,ivperp,iz]
        counter += 1
    end

    @begin_anyzv_region()
    @anyzv_serial_region begin
        @timeit_debug global_timer "ldiv!" ldiv!(this_output_buffer, precon_lu, this_input_buffer)
    end

    @begin_anyzv_region()
    counter = 1
    @loop_z_vperp_vpa iz ivperp ivpa begin
        precon_f[ivpa,ivperp,iz] = this_output_buffer[counter]
        counter += 1
    end

    # Ensure values of precon_f and precon_p are consistent across distributed-MPI block
    # boundaries. For precon_f take the upwind value, and for precon_p take the average.
    f_lower_endpoints = @view scratch_dummy.buffer_vpavperprs_1[:,:,ir,is]
    f_upper_endpoints = @view scratch_dummy.buffer_vpavperprs_2[:,:,ir,is]
    receive_buffer1 = @view scratch_dummy.buffer_vpavperprs_3[:,:,ir,is]
    receive_buffer2 = @view scratch_dummy.buffer_vpavperprs_4[:,:,ir,is]
    @begin_anyzv_vperp_vpa_region()
    @loop_vperp_vpa ivperp ivpa begin
        f_lower_endpoints[ivpa,ivperp] = precon_f[ivpa,ivperp,1]
        f_upper_endpoints[ivpa,ivperp] = precon_f[ivpa,ivperp,end]
    end
    # We upwind the z-derivatives in `z_advection_no_rs!()`, so would expect that
    # upwinding the results here in z would make sense.
    # However, (the following was true for the kinetic electron solve, but have not
    # checked again for the ion parallel dynamics solve) upwinding here makes convergence
    # much slower (~10x), compared to picking the values from one side or other of the
    # block boundary, or taking the average of the values on either side.  Neither
    # direction is special, so taking the average seems most sensible (although in an
    # intial test it does not seem to converge faster than just picking one or the other).
    # Maybe this could indicate that it is more important to have a fully self-consistent
    # Jacobian inversion for the `vpa_advection_no_rs()` part rather than taking half(ish)
    # of the values from one block and the other half(ish) from the other.
    reconcile_element_boundaries_MPI_z_pdf_vpavperpz!(
        precon_f, f_lower_endpoints, f_upper_endpoints, receive_buffer1,
        receive_buffer2, z)

    return nothing
end

function get_ion_parallel_dynamics_lu_preconditioner(
             nl_solver_params, f_new, buffer_1, buffer_2, buffer_3, buffer_4, dens, upar,
             p, dens_neutral, uz_neutral, p_neutral, moments, fields, composition,
             collisions, geometry, z::coordinate, vperp::coordinate, vpa::coordinate,
             z_spectral, vperp_spectral, vpa_spectral, r_advect, alpha_advect, z_advect,
             vperp_advect, vpa_advect, scratch_dummy, external_source_settings,
             num_diss_params, t_params, is::mk_int, ir::mk_int)

    recalculate_lu_preconditioner! = ion_parallel_dynamics_recalculate_lu!(
                                         nl_solver_params, f_new, dens, upar, p,
                                         dens_neutral, uz_neutral, p_neutral, moments,
                                         fields, composition, collisions, geometry, z,
                                         vperp, vpa, z_spectral, vperp_spectral,
                                         vpa_spectral, r_advect, alpha_advect, z_advect,
                                         vperp_advect, vpa_advect, scratch_dummy,
                                         external_source_settings, num_diss_params,
                                         t_params, is, ir)

    lu_precon! = ion_parallel_dynamics_lu_precon!(; nl_solver_params, scratch_dummy,
                                                  buffer_1, buffer_2, buffer_3, buffer_4,
                                                  z, is, ir)

    left_preconditioner = identity
    right_preconditioner = lu_precon!

    return left_preconditioner, right_preconditioner, recalculate_lu_preconditioner!
end

function get_ion_parallel_dynamics_preconditioner(
             nl_solver_params, f_new, buffer_1, buffer_2, buffer_3, buffer_4, dens, upar,
             p, dens_neutral, uz_neutral, p_neutral, moments, fields, composition,
             collisions, geometry, z, vperp, vpa, z_spectral, vperp_spectral,
             vpa_spectral, r_advect, alpha_advect, z_advect, vperp_advect, vpa_advect,
             scratch_dummy, external_source_settings, num_diss_params, t_params, is, ir)
    if nl_solver_params.preconditioner_type === Val(:ion_parallel_dynamics_lu)
        return get_ion_parallel_dynamics_lu_preconditioner(
                   nl_solver_params, f_new, buffer_1, buffer_2, buffer_3, buffer_4, dens,
                   upar, p, dens_neutral, uz_neutral, p_neutral, moments, fields,
                   composition, collisions, geometry, z, vperp, vpa, z_spectral,
                   vperp_spectral, vpa_spectral, r_advect, alpha_advect, z_advect,
                   vperp_advect, vpa_advect, scratch_dummy, external_source_settings,
                   num_diss_params, t_params, is, ir)
    elseif nl_solver_params.preconditioner_type === Val(:none)
        left_preconditioner = identity
        right_preconditioner = identity
        # Nothing to recalculate.
        recalculate_none_preconditioner = ()->nothing

        return left_preconditioner, right_preconditioner, recalculate_none_preconditioner
    else
        error("preconditioner_type=$(nl_solver_params.preconditioner_type) is not "
              * "supported for implicit ion parallel dynamics.")
    end
end

@kwdef @concrete struct ion_parallel_dynamics_residual!
    fvec_old
    pdf
    moments
    fields
    composition
    collisions
    geometry
    external_source_settings
    manufactured_source_list
    fp_arrays
    t_params
    scratch_dummy
    r
    z
    vperp
    vpa
    vzeta
    vr
    vz
    gyrophase
    spectral_objects
    advection_structs
    num_diss_params
    advance_implicit
    solve_counter
    dt::mk_float
    is::mk_int
    ir::mk_int
end

function (res::ion_parallel_dynamics_residual!)(this_residual, f_new; krylov=false)
    fvec_old = res.fvec_old
    pdf = res.pdf
    moments = res.moments
    fields = res.fields
    composition = res.composition
    collisions = res.collisions
    geometry = res.geometry
    external_source_settings = res.external_source_settings
    manufactured_source_list = res.manufactured_source_list
    fp_arrays = res.fp_arrays
    t_params = res.t_params
    scratch_dummy = res.scratch_dummy
    r = res.r
    z = res.z
    vperp = res.vperp
    vpa = res.vpa
    vzeta = res.vzeta
    vr = res.vr
    vz = res.vz
    gyrophase = res.gyrophase
    spectral_objects = res.spectral_objects
    advection_structs = res.advection_structs
    num_diss_params = res.num_diss_params
    advance_implicit = res.advance_implicit
    solve_counter = res.solve_counter
    dt = res.dt
    is = res.is
    ir = res.ir

    f_old = fvec_old.pdf
    density = fvec_old.density
    upar = fvec_old.upar
    p = fvec_old.p
    f_neutral = fvec_old.pdf_neutral
    density_neutral = fvec_old.density_neutral
    uz_neutral = fvec_old.uz_neutral
    p_neutral = fvec_old.p_neutral
    vperp_spectral = spectral_objects.vperp_spectral
    vpa_spectral = spectral_objects.vpa_spectral
    z_advect = @view advection_structs.z_advect[:,:,:,ir,is]
    vperp_advect = @view advection_structs.vperp_advect[:,:,:,ir,is]
    vpa_advect = @view advection_structs.vpa_advect[:,:,:,ir,is]

    enforce_boundary_conditions_f_no_sr!(
        f_new, density, upar, p, fields, moments, vpa, vperp, z, vpa_spectral,
        vperp_spectral, vpa_advect, vperp_advect, z_advect, composition, geometry,
        scratch_dummy, advance_implicit.vpa_diffusion, advance_implicit.vperp_diffusion,
        is, ir)

    @begin_anyzv_z_vperp_vpa_region()
    @loop_z_vperp_vpa iz ivperp ivpa begin
        this_residual[ivpa,ivperp,iz] = f_old[ivpa,ivperp,iz]
    end

    this_fvec_out = (pdf=this_residual, density=density, upar=upar, p=p,
                     pdf_neutral=f_neutral, density_neutral=density_neutral,
                     uz_neutral=uz_neutral, p_neutral=p_neutral)
    this_fvec_in = (pdf=f_new, density=density, upar=upar, p=p, pdf_neutral=f_neutral,
                    density_neutral=density_neutral, uz_neutral=uz_neutral,
                    p_neutral=p_neutral)
    ion_parallel_dynamics_euler_time_advance_no_sr!(
        this_fvec_out, this_fvec_in, pdf, nothing, fields, moments,
        advection_structs, vz, vr, vzeta, vpa, vperp, gyrophase, z, r, t_params,
        dt, spectral_objects, composition, collisions, geometry, scratch_dummy,
        manufactured_source_list, external_source_settings, num_diss_params,
        advance_implicit, fp_arrays, solve_counter[], is, ir)

    # Now
    #   residual = f_electron_old + dt*RHS(f_electron_newvar)
    # so update to desired residual
    @begin_anyzv_z_vperp_vpa_region()
    @loop_z_vperp_vpa iz ivperp ivpa begin
        this_residual[ivpa,ivperp,iz] = f_new[ivpa,ivperp,iz] - this_residual[ivpa,ivperp,iz]
    end

    # Set residual to zero where pdf_electron is determined by boundary conditions.
    if vpa.n > 1
        @begin_anyzv_z_vperp_region()
        @loop_z_vperp iz ivperp begin
            @views enforce_v_boundary_condition_local!(this_residual[:,ivperp,iz], vpa.bc,
                                                       vpa_advect[:,ivperp,iz],
                                                       advance_implicit.vpa_diffusion,
                                                       vpa, vpa_spectral)
        end
    end
    if vperp.n > 1
        @begin_anyzv_z_vpa_region()
        enforce_vperp_boundary_condition!(this_residual, vperp.bc, vperp, vperp_spectral,
                                          vperp_adv, advance_implicit.vperp_diffusion, ir)
    end
    zero_ion_parallel_dynamics_z_boundary_condition_points(this_residual, z, vpa, moments,
                                                           is, ir)

    return nothing
end

end # implicit_ion_parallel_dynamics
