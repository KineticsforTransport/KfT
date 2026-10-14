module ion_preconditioners

export get_ion_preconditioners

function get_ion_preconditioners(preconditioner_type, coords, outer_coords, spectral;
                                 boundary_skip_funcs::BSF=nothing) where {BSF}
    coord_sizes = Tuple(isa(c, coordinate) ? c.n : c for c ∈ coords)
    total_size_coords = prod(coord_sizes)
    outer_coord_sizes = Tuple(isa(c, coordinate) ? c.n : c for c ∈ outer_coords)

    if preconditioner_type === Val(:ion_parallel_lu)
        # One constraint for third_moment that feeds pdf into qpar.
        pdf_plus_constraints_size = total_size_coords + coords.z.n
        if coords.vperp.n == 1
            preconditioners = [(lu(sparse(1.0*I, 1, 1)),
                                create_jacobian_info(coords, spectral;
                                                     comm=comm_anyzv_subblock[],
                                                     synchronize=_anyzv_subblock_synchronize,
                                                     boundary_skip_funcs=boundary_skip_funcs.full,
                                                     ion_pdf=((:anyzv,:z,:vperp,:vpa), (:vpa, :vperp, :z), false),
                                                     third_moment=((:anyzv,:z), (:z,), true)),
                                allocate_shared_float(:newton_size=>pdf_plus_constraints_size;
                                                      comm=comm_anyzv_subblock[]),
                                allocate_shared_float(:newton_size=>pdf_plus_constraints_size;
                                                      comm=comm_anyzv_subblock[]),
                               )
                               for _ ∈ CartesianIndices(reverse(outer_coord_sizes))]
        else
            # 2V run, need another constraint for wpa^2 moment that feeds pdf into ppar.
            pdf_plus_constraints_size += coords.z.n

            preconditioners = [(lu(sparse(1.0*I, 1, 1)),
                                create_jacobian_info(coords, spectral;
                                                     comm=comm_anyzv_subblock[],
                                                     synchronize=_anyzv_subblock_synchronize,
                                                     boundary_skip_funcs=boundary_skip_funcs.full,
                                                     ion_pdf=((:anyzv,:z,:vperp,:vpa), (:vpa, :vperp, :z), false),
                                                     wpa2_moment=((:anyzv,:z), (:z,), true),
                                                     third_moment=((:anyzv,:z), (:z,), true)),
                                allocate_shared_float(:newton_size=>pdf_plus_constraints_size;
                                                      comm=comm_anyzv_subblock[]),
                                allocate_shared_float(:newton_size=>pdf_plus_constraints_size;
                                                      comm=comm_anyzv_subblock[]),
                               )
                               for _ ∈ CartesianIndices(reverse(outer_coord_sizes))]
        end
        # Initialise input buffers to zero so that constraint equations have zero on the
        # RHS.
        @begin_serial_region()
        @serial_region begin
            for p ∈ preconditioners
                p[3] .= 0.0
            end
        end
    elseif preconditioner_type === Val(:none)
        preconditioners = nothing
    else
        error("Unrecognised preconditioner_type=$preconditioner_type")
    end

    return preconditioners
end

@kwdef @concrete struct ion_parallel_dynamics_recalculate_lu!
    nl_solver_params
    f_new
    dens
    upar
    p
    dens_neutral
    moments
    fields
    collisions
    z
    vperp
    vpa
    z_spectral
    vperp_spectral
    vpa_spectral
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
    moments = recalc.moments
    fields = recalc.fields
    collisions = recalc.collisions
    z = recalc.z
    vperp = recalc.vperp
    vpa = recalc.vpa
    z_spectral = recalc.z_spectral
    vperp_spectral = recalc.vperp_spectral
    vpa_spectral = recalc.vpa_spectral
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

    nl_solver_params.solves_since_precon_update[] = 0
    nl_solver_params.precon_dt[] = t_params.dt[]

    orig_lu, precon, input_buffer, output_buffer =
        nl_solver_params.preconditioners[ir]

    fill_ion_kinetic_equation_Jacobian_evolve_nup!(
        precon, f_new, dens, upar, p, dens_neutral, moments, fields, collisions,
        external_source_settings, geometry, num_diss_params, scratch_dummy, r_advect,
        alpha_advect, z_advect, vperp_advect, vpa_advect, z, vperp, vpa, z_spectral,
        vperp_spectral, vpa_spectral, dt, is, ir)

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

@timeit_debug global_timer (lu_pre::ion_parallel_dynamics_lu_precon!)(x) = begin
    nl_solver_params = lu_pre.nl_solver_params
    scratch_dummy = lu_pre.scratch_dummy
    buffer_1 = lu_pre.buffer_1
    buffer_2 = lu_pre.buffer_2
    buffer_3 = lu_pre.buffer_3
    buffer_4 = lu_pre.buffer_4
    z = lu_pre.z
    is = lu_pre.is
    ir = lu_pre.ir

    precon_p, precon_f = x

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
             p, dens_neutral, moments, fields, collisions, z::coordinate,
             vperp::coordinate, vpa::coordinate, z_spectral, vperp_spectral, vpa_spectral,
             z_advect, vperp_advect, vpa_advect, scratch_dummy, external_source_settings,
             num_diss_params, t_params, is::mk_int, ir::mk_int)

    recalculate_lu_preconditioner! = ion_parallel_dynamics_recalculate_lu!(
                                         nl_solver_params, f_new, dens, upar, p,
                                         dens_neutral, moments, fields, collisions, z,
                                         vperp, vpa, z_spectral, vperp_spectral,
                                         vpa_spectral, z_advect, vperp_advect, vpa_advect,
                                         scratch_dummy, external_source_settings,
                                         num_diss_params, t_params, is, ir)

    lu_precon! = ion_parallel_dynamics_lu_precon!(; nl_solver_params, scratch_dummy,
                                                  buffer_1, buffer_2, buffer_3, buffer_4,
                                                  z, is, ir)

    left_preconditioner = identity
    right_preconditioner = lu_precon!

    return left_preconditioner, right_preconditioner, recalculate_lu_preconditioner!
end

function get_ion_parallel_dynamics_preconditioner(
             nl_solver_params, f_new, dens, upar, p, dens_neutral, moments, fields,
             collisions, z, vperp, vpa, z_spectral, vperp_spectral, vpa_spectral,
             z_advect, vperp_advect, vpa_advect, scratch_dummy, external_source_settings,
             num_diss_params, t_params, is, ir)
    if nl_solver_params.preconditioner_type === Val(:ion_parallel_dynamics_lu)
        return get_ion_parallel_dynamics_lu_preconditioner(
                   nl_solver_params, f_new, dens, upar, p, dens_neutral, moments, fields,
                   collisions, z, vperp, vpa, z_spectral, vperp_spectral, vpa_spectral,
                   z_advect, vperp_advect, vpa_advect, scratch_dummy,
                   external_source_settings, num_diss_params, t_params, is, ir)
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

end # ion_preconditioners
