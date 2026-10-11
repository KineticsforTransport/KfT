"""
"""
module source_terms

export source_terms!
export source_terms_manufactured!

using ..calculus: derivative!
using ..looping
using ..moment_kinetics_structs
using ..timer_utils

"""
calculate the source terms due to redefinition of the pdf to split off density,
flow and/or pressure, and use them to update the pdf
"""
@timeit global_timer source_terms!(
                         pdf_out, fvec_in, moments, r_advect, alpha_advect, z_advect, vpa,
                         vperp, z, r, dt, composition, collisions,
                         ion_source_settings) = begin
    evolve_density = moments.evolve_density
    evolve_upar = moments.evolve_upar
    evolve_p = moments.evolve_p
    if !evolve_density && !evolve_upar && !evolve_p
        return nothing
    end
    return source_terms!(pdf_out, fvec_in, moments, r_advect, alpha_advect, z_advect, vpa,
                         vperp, z, r, dt, composition, collisions, ion_source_settings,
                         Val(evolve_p))
end
function source_terms!(pdf_out, fvec_in, moments, r_advect, alpha_advect, z_advect, vpa,
                       vperp, z, r, dt, composition, collisions, ion_source_settings,
                       evolve_p::Val)

    @begin_s_r_z_vperp_vpa_region()

    pdf_in = fvec_in.pdf
    inner_args = get_source_terms_inner_args(vperp, moments.ion.vth, moments.ion.dvth_dt,
                                             moments.ion.dvth_dr, moments.ion.dvth_dz,
                                             fvec_in.density, moments.ion.ddens_dt,
                                             moments.ion.ddens_dr, moments.ion.ddens_dz,
                                             r_advect, alpha_advect, z_advect, pdf_in, dt,
                                             evolve_p)
    @loop_s_r is ir begin
        inner_args_sr = get_source_terms_inner_args_sr(is, ir, inner_args...)
        @loop_z iz begin
            inner_args_z = get_source_terms_inner_args_z(iz, inner_args_sr...)
            @loop_vperp_vpa ivperp ivpa begin
                pdf_out[ivpa,ivperp,iz,ir,is] += get_source_terms_inner(ivperp, ivpa,
                                                                        inner_args_z...)
            end
        end
    end

    return nothing
end

function source_terms_no_sr!(pdf_out, fvec_in, moments, r_advect, alpha_advect, z_advect,
                             vpa, vperp, z, r, dt, composition, collisions,
                             ion_source_settings, is, ir)
    evolve_density = moments.evolve_density
    evolve_upar = moments.evolve_upar
    evolve_p = moments.evolve_p
    if !evolve_density && !evolve_upar && !evolve_p
        return nothing
    end
    return source_terms_no_sr!(pdf_out, fvec_in, moments, r_advect, alpha_advect,
                               z_advect, vpa, vperp, z, r, dt, composition, collisions,
                               ion_source_settings, Val(evolve_p), is, ir)
end
function source_terms_no_sr!(pdf_out, fvec_in, moments, r_advect, alpha_advect, z_advect,
                             vpa, vperp, z, r, dt, composition, collisions,
                             ion_source_settings, evolve_p::Val, is, ir)

    @begin_anyzv_z_vperp_vpa_region()

    pdf_in = fvec_in.pdf
    @views inner_args =
        get_source_terms_inner_args(vperp, moments.ion.vth[:,ir,is],
                                    moments.ion.dvth_dt[:,ir,is],
                                    moments.ion.dvth_dr[:,ir,is],
                                    moments.ion.dvth_dz[:,ir,is], fvec_in.density,
                                    moments.ion.ddens_dt[:,ir,is],
                                    moments.ion.ddens_dr[:,ir,is],
                                    moments.ion.ddens_dz[:,ir,is], r_advect[:,:,:,ir,is],
                                    alpha_advect[:,:,:,ir,is], z_advect[:,:,:,ir,is],
                                    pdf_in, dt, evolve_p)
    @loop_z iz begin
        inner_args_z = get_source_terms_inner_args_z(iz, inner_args...)
        @loop_vperp_vpa ivperp ivpa begin
            pdf_out[ivpa,ivperp,iz] += get_source_terms_inner(ivperp, ivpa,
                                                              inner_args_z...)
        end
    end

    return nothing
end

@inline function get_source_terms_inner_args(vperp, vth, dvth_dt, dvth_dr, dvth_dz, n,
                                             dn_dt, dn_dr, dn_dz, r_advect, alpha_advect,
                                             z_advect, pdf_in, dt, evolve_p::Val{true})
    if vperp.n == 1
        # velocity dimension coefficient is 1.0 for 1V, and 3.0 for 3V. Only needed for evolve_p
        # since vth does not come in to gdot (Fdot) for other cases.
        v_dim_coeff = 1.0
    else
        v_dim_coeff = 3.0
    end
    return v_dim_coeff, vth, dvth_dt, dvth_dr, dvth_dz, n, dn_dt, dn_dr, dn_dz, r_advect,
           alpha_advect, z_advect, pdf_in, dt, evolve_p
end
@inline function get_source_terms_inner_args(vperp, vth, dvth_dt, dvth_dr, dvth_dz, n,
                                             dn_dt, dn_dr, dn_dz, r_advect, alpha_advect,
                                             z_advect, pdf_in, dt, evolve_p::Val{false})
    return n, dn_dt, dn_dr, dn_dz, r_advect, alpha_advect, z_advect, pdf_in, dt, evolve_p
end

@inline function get_source_terms_inner_args_sr(is, ir, v_dim_coeff, vth, dvth_dt,
                                                dvth_dr, dvth_dz, n, dn_dt, dn_dr, dn_dz,
                                                r_advect, alpha_advect, z_advect, pdf_in,
                                                dt, evolve_p::Val{true})
    return @views v_dim_coeff, vth[:,ir,is], dvth_dt[:,ir,is], dvth_dr[:,ir,is],
                  dvth_dz[:,ir,is], n[:,ir,is], dn_dt[:,ir,is], dn_dr[:,ir,is],
                  dn_dz[:,ir,is], r_advect[:,:,:,ir,is], alpha_advect[:,:,:,ir,is],
                  z_advect[:,:,:,ir,is], pdf_in[:,:,:,ir,is], dt, evolve_p
end
@inline function get_source_terms_inner_args_sr(is, ir, n, dn_dt, dn_dr, dn_dz, r_advect,
                                                alpha_advect, z_advect, pdf_in, dt,
                                                evolve_p::Val{false})
    return @views n[:,ir,is], dn_dt[:,ir,is], dn_dr[:,ir,is], dn_dz[:,ir,is],
                  r_advect[:,:,:,ir,is], alpha_advect[:,:,:,ir,is], z_advect[:,:,:,ir,is],
                  pdf_in[:,:,:,ir,is], dt, evolve_p
end

@inline function get_source_terms_inner_args_z(iz, v_dim_coeff, vth, dvth_dt, dvth_dr,
                                               dvth_dz, n, dn_dt, dn_dr, dn_dz, r_advect,
                                               alpha_advect, z_advect, pdf_in, dt,
                                               evolve_p::Val{true})
    ddt_term = v_dim_coeff / vth[iz] * dvth_dt[iz] - dn_dt[iz] / n[iz]
    rdot_coefficient = v_dim_coeff / vth[iz] * dvth_dr[iz] - dn_dr[iz] / n[iz]
    zdot_coefficient = v_dim_coeff / vth[iz] * dvth_dz[iz] - dn_dz[iz] / n[iz]
    return @views ddt_term, rdot_coefficient, zdot_coefficient, r_advect[:,:,iz],
                  alpha_advect[:,:,iz], z_advect[:,:,iz], pdf_in[:,:,iz], dt
end
@inline function get_source_terms_inner_args_z(iz, n, dn_dt, dn_dr, dn_dz, r_advect,
                                               alpha_advect, z_advect, pdf_in, dt,
                                               evolve_p::Val{false})
    ddt_term = - dn_dt[iz] / n[iz]
    rdot_coefficient = - dn_dr[iz] / n[iz]
    zdot_coefficient = - dn_dz[iz] / n[iz]
    return @views ddt_term, rdot_coefficient, zdot_coefficient, r_advect[:,:,iz],
                  alpha_advect[:,:,iz], z_advect[:,:,iz], pdf_in[:,:,iz], dt
end

@inline function get_source_terms_inner(ivperp, ivpa, ddt_term, rdot_coefficient,
                                        zdot_coefficient, r_advect, alpha_advect,
                                        z_advect, pdf_in, dt)
    return dt * (
                 ddt_term
                 + rdot_coefficient * r_advect[ivpa,ivperp]
                 + zdot_coefficient * (alpha_advect[ivpa,ivperp] + z_advect[ivpa,ivperp])
                ) *
           pdf_in[ivpa,ivperp]
end

function get_source_terms_evolve_nup(sub_terms::IonSubTerms)
    nvperp = sub_terms.nvperp
    vth = sub_terms.vth
    dvth_dt = sub_terms.dvth_dt
    dvth_dr = sub_terms.dvth_dr
    dvth_dz = sub_terms.dvth_dz
    n = sub_terms.n
    dn_dt = sub_terms.dn_dt
    dn_dr = sub_terms.dn_dr
    dn_dz = sub_terms.dn_dz
    r_speed = sub_terms.r_speed
    alpha_speed = sub_terms.alpha_speed
    z_speed = sub_terms.z_speed
    f = sub_terms.f

    # velocity dimension coefficient is 1.0 for 1V, and 3.0 for 3V. Only needed for
    # evolve_p since vth does not come in to gdot (Fdot) for other cases.
    if nvperp == 1
        v_dim_coeff = 1.0
    else
        v_dim_coeff = 3.0
    end

    ddt_term = v_dim_coeff * vth^(-1) * dvth_dt - dn_dt * n^(-1)
    rdot_coefficient = v_dim_coeff * vth^(-1) * dvth_dr - dn_dr * n^(-1)
    zdot_coefficient = v_dim_coeff * vth^(-1) * dvth_dz - dn_dz * n^(-1)
    term = -(ddt_term
             + rdot_coefficient * r_speed
             + zdot_coefficient * (alpha_speed + z_speed)
            ) * f

    return term
end

"""
calculate the source terms due to redefinition of the pdf to split off density,
flow and/or pressure, and use them to update the pdf
"""
@timeit global_timer source_terms_neutral!(
                         pdf_out, fvec_in, moments, vz, z, r, dt, composition, collisions,
                         neutral_source_settings) = begin
    evolve_density = moments.evolve_density
    evolve_upar = moments.evolve_upar
    evolve_p = moments.evolve_p
    if !evolve_density && !evolve_upar && !evolve_p
        return nothing
    end
    return source_terms_neutral!(pdf_out, fvec_in, moments, vz, z, r, dt, composition,
                                 collisions, neutral_source_settings, Val(evolve_upar),
                                 Val(evolve_p))
end
function source_terms_neutral!(pdf_out, fvec_in, moments, vz, z, r, dt, composition,
                               collisions, neutral_source_settings, evolve_upar::Val,
                               evolve_p::Val)

    @begin_sn_r_z_vzeta_vr_vz_region()

    inner_args =
        get_source_terms_neutral_inner_args(vz, moments.neutral.vth,
                                            moments.neutral.dvth_dt,
                                            moments.neutral.dvth_dz, fvec_in.uz_neutral,
                                            fvec_in.density_neutral,
                                            moments.neutral.ddens_dt,
                                            moments.neutral.ddens_dz, fvec_in.pdf_neutral,
                                            dt, evolve_upar, evolve_p)

    @loop_sn_r isn ir begin
        inner_args_snr = get_source_terms_neutral_inner_args_snr(isn, ir, inner_args...)
        @loop_z iz begin
            inner_args_z = get_source_terms_neutral_inner_args_z(iz, inner_args_snr...)
            @loop_vzeta_vr_vz ivzeta ivr ivz begin
                pdf_out[ivz,ivr,ivzeta,iz,ir,isn] +=
                    get_source_terms_neutral_inner(ivzeta, ivr, ivz, inner_args_z...)
            end
        end
    end

    return nothing
end

@inline function get_source_terms_neutral_inner_args(vz, vth, dvth_dt, dvth_dz, uz, n,
                                                     dn_dt, dn_dz, pdf_in, dt,
                                                     evolve_upar::Val{true},
                                                     evolve_p::Val{true})
    return vz.grid, vth, dvth_dt, dvth_dz, uz, n, dn_dt, dn_dz, pdf_in, dt, evolve_upar,
           evolve_p
end
@inline function get_source_terms_neutral_inner_args(vz, vth, dvth_dt, dvth_dz, uz, n,
                                                     dn_dt, dn_dz, pdf_in, dt,
                                                     evolve_upar::Val{true},
                                                     evolve_p::Val{false})
    return vz.grid, uz, n, dn_dt, dn_dz, pdf_in, dt, evolve_upar, evolve_p
end
@inline function get_source_terms_neutral_inner_args(vz, vth, dvth_dt, dvth_dz, uz, n,
                                                     dn_dt, dn_dz, pdf_in, dt,
                                                     evolve_upar::Val{false},
                                                     evolve_p::Val{false})
    return vz.grid, n, dn_dt, dn_dz, pdf_in, dt, evolve_upar, evolve_p
end

@inline function get_source_terms_neutral_inner_args_snr(isn, ir, vz_grid, vth, dvth_dt,
                                                         dvth_dz, uz, n, dn_dt, dn_dz,
                                                         pdf_in, dt,
                                                         evolve_upar::Val{true},
                                                         evolve_p::Val{true})
    return @views vz_grid, vth[:,ir,isn], dvth_dt[:,ir,isn], dvth_dz[:,ir,isn],
                  uz[:,ir,isn], n[:,ir,isn], dn_dt[:,ir,isn], dn_dz[:,ir,isn],
                  pdf_in[:,:,:,:,ir,isn], dt, evolve_upar, evolve_p
end
@inline function get_source_terms_neutral_inner_args_snr(isn, ir, vz_grid, uz, n, dn_dt, dn_dz,
                                                         pdf_in, dt,
                                                         evolve_upar::Val{true},
                                                         evolve_p::Val{false})
    return @views vz_grid, uz[:,ir,isn], n[:,ir,isn], dn_dt[:,ir,isn],
                  dn_dz[:,ir,isn], pdf_in[:,:,:,:,ir,isn], dt, evolve_upar, evolve_p
end
@inline function get_source_terms_neutral_inner_args_snr(isn, ir, vz_grid, n, dn_dt, dn_dz,
                                                         pdf_in, dt,
                                                         evolve_upar::Val{false},
                                                         evolve_p::Val{false})
    return @views vz_grid, n[:,ir,isn], dn_dt[:,ir,isn], dn_dz[:,ir,isn],
                  pdf_in[:,:,:,:,ir,isn], dt, evolve_upar, evolve_p
end

@inline function get_source_terms_neutral_inner_args_z(iz, vz_grid, vth, dvth_dt, dvth_dz, uz,
                                                       n, dn_dt, dn_dz, pdf_in, dt,
                                                       evolve_upar::Val{true},
                                                       evolve_p::Val{true})
    coefficient1 = -(dn_dt[iz] + uz[iz] * dn_dz[iz]) / n[iz] +
                   (dvth_dt[iz] + uz[iz] * dvth_dz[iz]) / vth[iz]
    coefficient2 = -vth[iz] * dn_dz[iz] / n[iz] + dvth_dz[iz]
    return @views vz_grid, coefficient1, coefficient2, pdf_in[:,:,:,iz], dt
end
@inline function get_source_terms_neutral_inner_args_z(iz, vz_grid, uz, n, dn_dt, dn_dz,
                                                       pdf_in, dt, evolve_upar::Val{true},
                                                       evolve_p::Val{false})
    coefficient1 = -(dn_dt[iz] + uz[iz] * dn_dz[iz]) / n[iz]
    coefficient2 = -dn_dz[iz] / n[iz]
    return @views vz_grid, coefficient1, coefficient2, pdf_in[:,:,:,iz], dt
end
@inline function get_source_terms_neutral_inner_args_z(iz, vz_grid, n, dn_dt, dn_dz,
                                                       pdf_in, dt,
                                                       evolve_upar::Val{false},
                                                       evolve_p::Val{false})
    coefficient1 = -dn_dt[iz] / n[iz]
    coefficient2 = -dn_dz[iz] / n[iz]
    return @views vz_grid, coefficient1, coefficient2, pdf_in[:,:,:,iz], dt
end

@inline function get_source_terms_neutral_inner(ivzeta, ivr, ivz, vz_grid, coefficient1,
                                                coefficient2, pdf_in, dt)
    return dt * (coefficient1 + vz_grid[ivz] * coefficient2) * pdf_in[ivz,ivr,ivzeta]
end

"""
advance the dfn with an arbitrary source function 
"""
@timeit global_timer source_terms_manufactured!(
                         pdf_ion_out, pdf_neutral_out, vz, vr, vzeta, vpa, vperp, z, r, t,
                         dt, composition, manufactured_source_list) = begin
    if manufactured_source_list.time_independent_sources
        # the (time-independent) manufactured source arrays
        Source_i = manufactured_source_list.Source_i_array
        Source_n = manufactured_source_list.Source_n_array

        @begin_s_r_z_region()

        @loop_s is begin
            @loop_r_z_vperp_vpa ir iz ivperp ivpa begin
                pdf_ion_out[ivpa,ivperp,iz,ir,is] += dt*Source_i[ivpa,ivperp,iz,ir]
            end
        end

        if composition.n_neutral_species > 0
            @begin_sn_r_z_region()
            @loop_sn isn begin
                @loop_r_z_vzeta_vr_vz ir iz ivzeta ivr ivz begin
                    pdf_neutral_out[ivz,ivr,ivzeta,iz,ir,isn] += dt*Source_n[ivz,ivr,ivzeta,iz,ir]
                end
            end
        end
    else
        # the manufactured source functions
        Source_i_func = manufactured_source_list.Source_i_func
        Source_n_func = manufactured_source_list.Source_n_func

        @begin_s_r_z_region()

        @loop_s is begin
            @loop_r_z_vperp_vpa ir iz ivperp ivpa begin
                pdf_ion_out[ivpa,ivperp,iz,ir,is] += dt*Source_i_func(vpa.grid[ivpa],vperp.grid[ivperp],z.grid[iz],r.grid[ir],t)
            end
        end

        if composition.n_neutral_species > 0
            @begin_sn_r_z_region()
            @loop_sn isn begin
                @loop_r_z_vzeta_vr_vz ir iz ivzeta ivr ivz begin
                    pdf_neutral_out[ivz,ivr,ivzeta,iz,ir,isn] += dt*Source_n_func(vz.grid[ivz],vr.grid[ivr],vzeta.grid[ivzeta],z.grid[iz],r.grid[ir],t)
                end
            end
        end
    end
    return nothing
end

end
