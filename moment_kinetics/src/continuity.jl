"""
"""
module continuity

export continuity_equation!

using ..calculus: derivative!
using ..looping
using ..timer_utils

"""
use the continuity equation dn/dt + d(n*upar)/dz to update the density n for all ion
species
"""
@timeit global_timer continuity_equation!(
                         dens_out, fvec_in, fields, moments, composition, geometry, dt,
                         ionization, ion_source_settings, num_diss_params) = begin
    @begin_s_r_z_region()

    ddens_dt = moments.ion.ddens_dt

    vEr = fields.vEr
    vEz = fields.vEz
    n = fvec_in.density
    upar = fvec_in.upar
    dn_dr_upwind = moments.ion.ddens_dr_upwind
    dn_dz_upwind = moments.ion.ddens_dz_upwind
    du_dz = moments.ion.dupar_dz
    bz = geometry.bzed
    Bmag = geometry.Bmag
    dBdz = geometry.dBdz
    dBdr = geometry.dBdr

    @loop_s_r_z is ir iz begin
        # ddens_dz is upwinded using upar
        ddens_dt[iz,ir,is] =
            get_ddens_dt_inner_main(n[iz,ir,is], upar[iz,ir,is], dn_dr_upwind[iz,ir,is],
                                    dn_dz_upwind[iz,ir,is], du_dz[iz,ir,is], vEr[iz,ir],
                                    vEz[iz,ir], bz[iz,ir], Bmag[iz,ir], dBdr[iz,ir],
                                    dBdz[iz,ir])
    end

    # update the density to account for ionization collisions;
    # ionization collisions increase the density for ions and decrease the density for neutrals
    if composition.n_neutral_species > 0 && ionization !== nothing
        density_neutral = fvec.density_neutral
        @loop_s_r_z is ir iz begin
            ddens_dt[iz,ir,is] +=
                get_ddens_dt_reactions_inner(ionization, density[iz,ir,is],
                                             density_neutral[iz,ir,is])
        end
    end

    for index ∈ eachindex(ion_source_settings)
        if ion_source_settings[index].active
            @views source_amplitude = moments.ion.external_source_density_amplitude[:, :, index]
            @loop_s_r_z is ir iz begin
                ddens_dt[iz,ir,is] += source_amplitude[iz,ir]
            end
        end
    end

    # Ad-hoc diffusion to stabilise numerics...
    diffusion_coefficient = num_diss_params.ion.moment_dissipation_coefficient
    if diffusion_coefficient > 0.0
        @loop_s_r_z is ir iz begin
            ddens_dt[iz,ir,is] += diffusion_coefficient*moments.ion.d2dens_dz2[iz,ir,is]
        end
    end

    @loop_s_r_z is ir iz begin
        dens_out[iz,ir,is] += dt * ddens_dt[iz,ir,is]
    end

    return nothing
end

# This version only calculates dnupar_dt, and does not update a 'new upar'.
function continuity_equation_no_sr!(
                         fvec_in, fields, moments, composition, geometry, dt,
                         ionization, ion_source_settings, num_diss_params, is, ir)
    @begin_anyzv_z_region()

    ddens_dt = @view moments.ion.ddens_dt[:,ir,is]

    vEr = @view fields.vEr[:,ir]
    vEz = @view fields.vEz[:,ir]
    n = @view fvec_in.density[:,ir,is]
    upar = @view fvec_in.upar[:,ir,is]
    dn_dr_upwind = @view moments.ion.ddens_dr_upwind[:,ir,is]
    dn_dz_upwind = @view moments.ion.ddens_dz_upwind[:,ir,is]
    du_dz = @view moments.ion.dupar_dz[:,ir,is]
    bz = @view geometry.bzed[:,ir]
    Bmag = @view geometry.Bmag[:,ir]
    dBdz = @view geometry.dBdz[:,ir]
    dBdr = @view geometry.dBdr[:,ir]

    @loop_z iz begin
        # ddens_dz is upwinded using upar
        ddens_dt[iz] =
            get_ddens_dt_inner_main(n[iz], upar[iz], dn_dr_upwind[iz], dn_dz_upwind[iz],
                                    du_dz[iz], vEr[iz], vEz[iz], bz[iz], Bmag[iz],
                                    dBdr[iz], dBdz[iz])
    end

    # update the density to account for ionization collisions;
    # ionization collisions increase the density for ions and decrease the density for neutrals
    if composition.n_neutral_species > 0 && ionization !== nothing
        density_neutral = @view fvec.density_neutral[:,ir,is]
        @loop_z iz begin
            ddens_dt[iz] +=
                get_ddens_dt_reactions_inner(ionization, density[iz], density_neutral[iz])
        end
    end

    for index ∈ eachindex(ion_source_settings)
        if ion_source_settings[index].active
            @views source_amplitude = moments.ion.external_source_density_amplitude[:,ir,index]
            @loop_z iz begin
                ddens_dt[iz] += source_amplitude[iz]
            end
        end
    end

    # Ad-hoc diffusion to stabilise numerics...
    diffusion_coefficient = num_diss_params.ion.moment_dissipation_coefficient
    if diffusion_coefficient > 0.0
        @loop_z iz begin
            ddens_dt[iz] += diffusion_coefficient*moments.ion.d2dens_dz2[iz]
        end
    end

    return nothing
end

@inline function get_ddens_dt_inner_main(n, upar, dn_dr_upwind, dn_dz_upwind, du_dz, vEr,
                                         vEz, bz, Bmag, dBdr, dBdz)
    return -(vEr * dn_dr_upwind
             + (vEz + bz * upar) * dn_dz_upwind
             + bz * n * du_dz
             - (1/Bmag) * (vEr * dBdr + (vEz + bz * upar) * dBdz) * n)
end

@inline function get_ddens_dt_reactions_inner(ionization, density, density_neutral)
    return ionization * density * density_neutral
end

"""
use the continuity equation dn/dt + d(n*upar)/dz to update the density n for all neutral
species
"""
@timeit global_timer neutral_continuity_equation!(
                         dens_out, fvec_in, moments, composition, dt, ionization,
                         neutral_source_settings, num_diss_params) = begin
    @begin_sn_r_z_region()

    ddens_dt = moments.neutral.ddens_dt

    @loop_sn_r_z isn ir iz begin
        # Use ddens_dz is upwinded using uz
        ddens_dt[iz,ir,isn] =
            - (fvec_in.uz_neutral[iz,ir,isn]*moments.neutral.ddens_dz_upwind[iz,ir,isn] +
               fvec_in.density_neutral[iz,ir,isn]*moments.neutral.duz_dz[iz,ir,isn])
    end

    # update the density to account for ionization collisions;
    # ionization collisions increase the density for ions and decrease the density for neutrals
    if composition.n_neutral_species > 0 && ionization !== nothing
        @loop_sn_r_z isn ir iz begin
            ddens_dt[iz,ir,isn] -= ionization*fvec_in.density[iz,ir,isn]*fvec_in.density_neutral[iz,ir,isn]
        end
    end

    for index ∈ eachindex(neutral_source_settings)
        if neutral_source_settings[index].active
            @views source_amplitude = moments.neutral.external_source_density_amplitude[:, :, index]
            @loop_s_r_z is ir iz begin
                ddens_dt[iz,ir,is] += source_amplitude[iz,ir]
            end
        end
    end

    # Ad-hoc diffusion to stabilise numerics...
    diffusion_coefficient = num_diss_params.neutral.moment_dissipation_coefficient
    if diffusion_coefficient > 0.0
        @loop_sn_r_z isn ir iz begin
            ddens_dt[iz,ir,isn] += diffusion_coefficient*moments.neutral.d2dens_dz2[iz,ir,isn]
        end
    end

    @loop_sn_r_z isn ir iz begin
        # Use ddens_dz is upwinded using uz
        dens_out[iz,ir,isn] += dt * ddens_dt[iz,ir,isn]
    end

    return nothing
end

end
