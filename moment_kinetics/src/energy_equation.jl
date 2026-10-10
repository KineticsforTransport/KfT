"""
"""
module energy_equation

export energy_equation!
export neutral_energy_equation!

using ..calculus: derivative!
using ..looping
using ..timer_utils

"""
evolve the parallel pressure by solving the energy equation
"""
@timeit global_timer energy_equation!(
                         p_out, fvec, moments, fields, collisions, dt, composition,
                         geometry, ion_source_settings, num_diss_params) = begin

    @begin_s_r_z_region()

    density = fvec.density
    upar = fvec.upar
    p = fvec.p
    ppar = moments.ion.ppar
    pperp = moments.ion.pperp
    dupar_dz = moments.ion.dupar_dz
    dp_dr_upwind = moments.ion.dp_dr_upwind
    dp_dz_upwind = moments.ion.dp_dz_upwind
    qpar = moments.ion.qpar
    dqpar_dz = moments.ion.dqpar_dz
    dp_dt = moments.ion.dp_dt
    vEr = fields.vEr
    vEz = fields.vEz
    bz = geometry.bzed
    Bmag = geometry.Bmag
    dBdz = geometry.dBdz
    dBdr = geometry.dBdr

    @loop_s_r_z is ir iz begin
        dp_dt[iz,ir,is] = get_dpdt_inner_main(upar[iz,ir,is], p[iz,ir,is], ppar[iz,ir,is],
                                              pperp[iz,ir,is], dupar_dz[iz,ir,is],
                                              dp_dr_upwind[iz,ir,is],
                                              dp_dz_upwind[iz,ir,is], dqpar_dz[iz,ir,is],
                                              bz[iz,ir], Bmag[iz,ir], dBdr[iz,ir],
                                              dBdz[iz,ir], vEr[iz,ir], vEz[iz,ir])
    end


    for index ∈ eachindex(ion_source_settings)
        if ion_source_settings[index].active
            @views source_amplitude = moments.ion.external_source_pressure_amplitude[:, :, index]
            @loop_s_r_z is ir iz begin
                dp_dt[iz,ir,is] += source_amplitude[iz,ir]
            end
        end
    end

    diffusion_coefficient = num_diss_params.ion.moment_dissipation_coefficient
    if diffusion_coefficient > 0.0
        @loop_s_r_z is ir iz begin
            dp_dt[iz,ir,is] += diffusion_coefficient*moments.ion.d2p_dz2[iz,ir,is]
        end
    end

    # add in contributions due to charge exchange/ionization collisions
    if composition.n_neutral_species > 0
        charge_exchange = collisions.reactions.charge_exchange_frequency
        ionization = collisions.reactions.ionization_frequency
        density_neutral = fvec.density_neutral
        uz_neutral = fvec.uz_neutral
        p_neutral = fvec.p_neutral
        if charge_exchange !== nothing || ionization !== nothing
            @loop_s_r_z is ir iz begin
                dp_dt[iz,ir,is] +=
                    get_ion_reactions_inner(charge_exchange, ionization,
                                            density[iz,ir,is], upar[iz,ir,is],
                                            p[iz,ir,is], density_neutral[iz,ir,is],
                                            uz_neutral[iz,ir,is], p_neutral[iz,ir,is])
            end
        end
    end

    @loop_s_r_z is ir iz begin
        p_out[iz,ir,is] += dt * dp_dt[iz,ir,is]
    end

    return nothing
end

# This version only calculates dp_dt, and does not update a 'new p'.
function energy_equation_no_sr!(fvec, moments, fields, collisions, dt, composition,
                                geometry, ion_source_settings, num_diss_params, is, ir)

    @begin_anyzv_z_region()

    density = @view fvec.density[:,ir,is]
    upar = @view fvec.upar[:,ir,is]
    p = @view fvec.p[:,ir,is]
    ppar = @view moments.ion.ppar[:,ir,is]
    pperp = @view moments.ion.pperp[:,ir,is]
    dupar_dz = @view moments.ion.dupar_dz[:,ir,is]
    dp_dr_upwind = @view moments.ion.dp_dr_upwind[:,ir,is]
    dp_dz_upwind = @view moments.ion.dp_dz_upwind[:,ir,is]
    qpar = @view moments.ion.qpar[:,ir,is]
    dqpar_dz = @view moments.ion.dqpar_dz[:,ir,is]
    dp_dt = @view moments.ion.dp_dt[:,ir,is]
    vEr = @view fields.vEr[:,ir]
    vEz = @view fields.vEz[:,ir]
    bz = @view geometry.bzed[:,ir]
    Bmag = @view geometry.Bmag[:,ir]
    dBdz = @view geometry.dBdz[:,ir]
    dBdr = @view geometry.dBdr[:,ir]

    @loop_z iz begin
        dp_dt[iz] = get_dpdt_inner_main(upar[iz], p[iz], ppar[iz], pperp[iz],
                                        dupar_dz[iz], dp_dr_upwind[iz], dp_dz_upwind[iz],
                                        dqpar_dz[iz], bz[iz], Bmag[iz], dBdr[iz],
                                        dBdz[iz], vEr[iz], vEz[iz])
    end


    for index ∈ eachindex(ion_source_settings)
        if ion_source_settings[index].active
            @views source_amplitude = moments.ion.external_source_pressure_amplitude[:,ir,index]
            @loop_z iz begin
                dp_dt[iz] += source_amplitude[iz]
            end
        end
    end

    diffusion_coefficient = num_diss_params.ion.moment_dissipation_coefficient
    if diffusion_coefficient > 0.0
        @loop_z iz begin
            dp_dt[iz] += diffusion_coefficient*moments.ion.d2p_dz2[iz]
        end
    end

    # add in contributions due to charge exchange/ionization collisions
    if composition.n_neutral_species > 0
        charge_exchange = collisions.reactions.charge_exchange_frequency
        ionization = collisions.reactions.ionization_frequency
        density_neutral = fvec.density_neutral
        uz_neutral = fvec.uz_neutral
        p_neutral = fvec.p_neutral
        if charge_exchange !== nothing || ionization !== nothing
            @loop_z iz begin
                dp_dt[iz] +=
                    get_dpdt_reactions_inner(charge_exchange, ionization, density[iz],
                                             upar[iz], p[iz], density_neutral[iz],
                                             uz_neutral[iz], p_neutral[iz])
            end
        end
    end

    return nothing
end

@inline function get_dpdt_inner_main(upar, p, ppar, pperp, dupar_dz, dp_dr_upwind,
                                     dp_dz_upwind, dqpar_dz, bz, Bmag, dBdr, dBdz, vEr,
                                     vEz)
    return -(vEr * dp_dr_upwind
             + (vEz + bz * upar) * dp_dz_upwind
             + bz * p * dupar_dz
             + 2.0/3.0 * bz * dqpar_dz
             + 2.0/3.0 * bz * ppar * dupar_dz
             - 2.0/3.0 * (1/Bmag) * ((2 * pperp + 0.5 * ppar) *
                                     (vEr * dBdr[iz,ir] + (vEz + bz * upar) * dBdz)
                                     + bz * dBdz * qpar))
end

@inline function get_dpdt_reactions_inner(charge_exchange, ionization, density, upar, p,
                                          density_neutral, uz_neutral, p_neutral)
    if charge_exchange !== nothing
        result =
            - charge_exchange * (
                 density_neutral * p - density * p_neutral -
                 1.0/3.0 * density * density_neutral * (upar - uz_neutral)^2)
    else
        result = 0.0
    end
    if ionization !== nothing
        result +=
            ionization*density * (
                p_neutral + 1.0/3.0 * density_neutral * (upar - uz_neutral)^2)
    end
    return result
end

"""
evolve the neutral parallel pressure by solving the energy equation
"""
@timeit global_timer neutral_energy_equation!(
                         p_out, fvec, moments, collisions, dt, composition,
                         neutral_source_settings, num_diss_params) = begin

    @begin_sn_r_z_region()

    dp_dt = moments.neutral.dp_dt

    @loop_sn_r_z isn ir iz begin
        dp_dt[iz,ir,isn] = (-fvec.uz_neutral[iz,ir,isn]*moments.neutral.dp_dz_upwind[iz,ir,isn]
                            -fvec.p_neutral[iz,ir,isn]*moments.neutral.duz_dz[iz,ir,isn]
                            - 2.0/3.0*moments.neutral.dqz_dz[iz,ir,isn]
                            - 2.0/3.0*moments.neutral.pz[iz,ir,isn]*moments.neutral.duz_dz[iz,ir,isn])
    end

    for index ∈ eachindex(neutral_source_settings)
        if neutral_source_settings[index].active
            @views source_amplitude = moments.neutral.external_source_pressure_amplitude[:, :, index]
            @loop_s_r_z isn ir iz begin
                dp_dt[iz,ir,isn] += source_amplitude[iz,ir]
            end
        end
    end

    diffusion_coefficient = num_diss_params.neutral.moment_dissipation_coefficient
    if diffusion_coefficient > 0.0
        @loop_sn_r_z isn ir iz begin
            dp_dt[iz,ir,isn] += diffusion_coefficient*moments.neutral.d2pz_dz2[iz,ir,isn]
        end
    end

    # add in contributions due to charge exchange/ionization collisions
    if composition.n_neutral_species > 0
        charge_exchange = collisions.reactions.charge_exchange_frequency
        ionization = collisions.reactions.ionization_frequency
        if charge_exchange !== nothing
            @loop_sn_r_z isn ir iz begin
                dp_dt[iz,ir,isn] -=
                    charge_exchange*(
                        fvec.density[iz,ir,isn]*fvec.p_neutral[iz,ir,isn] -
                        fvec.density_neutral[iz,ir,isn]*fvec.p[iz,ir,isn] -
                        1.0/3.0 * fvec.density_neutral[iz,ir,isn]*fvec.density[iz,ir,isn] *
                            (fvec.uz_neutral[iz,ir,isn] - fvec.upar[iz,ir,isn])^2)
            end
        end
        if ionization !== nothing
            @loop_sn_r_z isn ir iz begin
                dp_dt[iz,ir,isn] -=
                    ionization*fvec.density[iz,ir,isn]*fvec.p_neutral[iz,ir,isn]
            end
        end
    end

    @loop_sn_r_z isn ir iz begin
        p_out[iz,ir,isn] += dt * dp_dt[iz,ir,isn]
    end

    return nothing
end

"""
    get_dvth_dt_expanded_term_evolve_nup(bzed, mi, n, vth, ppar, dupar_dz, dqpar_dz)

Note that this function only includes terms that depend directly on the ion shape
function, as these are the only ones that will contribute to the Jacobian matrix from
here. The contributions of all the 'constant' terms contribute only through the
`dupar_dt_array`.
"""
function get_dvth_dt_expanded_term_evolve_nup(bzed, mi, n, vth, ppar, dupar_dz, dqpar_dz)
    return - 2.0 * (3.0 * mi * n * vth)^(-1) * bzed * (dqpar_dz + ppar * dupar_dz)
end

end
