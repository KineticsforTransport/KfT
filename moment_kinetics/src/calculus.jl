"""
"""
module calculus

# Import moment_kinetics so that we can refer to it in docstrings
import moment_kinetics

export derivative!, second_derivative!, laplacian_derivative!
export elementwise_indefinite_integration!
export reconcile_element_boundaries_MPI!
export integral
export indefinite_integral!

using ..array_allocation: allocate_float
using ..moment_kinetics_structs: discretization_info, null_spatial_dimension_info,
                                 null_velocity_dimension_info, null_vperp_dimension_info,
                                 weak_discretization_info, fourier_info
using ..timer_utils
using ..type_definitions: mk_float, mk_int
using MPI
using ..communication: block_rank
using ..communication: @_block_synchronize, @_anyzv_subblock_synchronize
using ..debugging
using ..looping

using LinearAlgebra

"""
    elementwise_derivative!(coord, f, speed, spectral)
    elementwise_derivative!(coord, f, spectral)

Generic function for element-by-element derivatives

First signature, with `speed`, calculates an upwind derivative, the second signature
calculates a derivative without upwinding information.

Result is stored in coord.scratch_2d.
"""
function elementwise_derivative! end

function elementwise_derivative!(coord, f, speed,
                                 null::Union{null_spatial_dimension_info,null_velocity_dimension_info,null_vperp_dimension_info})
    coord.scratch_2d .= 0.0
    return nothing
end

"""
"""
function elementwise_indefinite_integration! end

"""
    indefinite_integral!(pf, f, coord, spectral)

Indefinite line integral. 

This function is designed to work on local-in-memory data only,
with distributed-memory MPI not implemented here.
A function which integrates along a line which is distributed
in memory exists in [`moment_kinetics.em_fields`](@ref) as 
`calculate_phi_from_Epar!()`. The distributed-memory functionality could
be ported to a generic function, similiar to how the derivative!
functions are generalised in [`moment_kinetics.derivatives`](@ref).
"""
function indefinite_integral!(pf, f, coord, spectral::discretization_info)
    # get the indefinite integral at each grid point within each element and store in
    # coord.scratch_2d
    elementwise_indefinite_integration!(coord, f, spectral)
    # map the integral from the elemental grid to the full grid;
    # taking care to match integral constants at the boundaries
    # we assume that the lower limit of the indefinite integral
    # is the lowest value in coord.grid
    indefinite_integral_elements_to_full_grid!(pf, coord)
end

"""
"""
function indefinite_integral_elements_to_full_grid!(pf, coord)
    pf2d = coord.scratch_2d
    nelement_local = coord.nelement_local
    ngrid = coord.ngrid
    igrid_full = coord.igrid_full   
    j = 1 # the first element
    for i in 1:ngrid
        pf[i] = pf2d[i,j]
    end
    for j in 2:nelement_local
        ilast = igrid_full[ngrid,j-1] # the grid index of the last point in the previous element
        for k in 2:coord.ngrid
            i = coord.igrid_full[k,j] # the index in the full grid
            pf[i] = pf2d[k,j] + pf[ilast]
        end
    end
    return nothing
end

"""
    derivative!(df, f, coord, speed, spectral)

Upwinding derivative.
"""
function derivative!(df, f, coord, speed, spectral::discretization_info)
    # get the derivative at each grid point within each element and store in
    # coord.scratch_2d
    elementwise_derivative!(coord, f, speed, spectral)
    # map the derivative from the elemental grid to the full grid;
    # at element boundaries, use the derivative from the upwind element.
    derivative_elements_to_full_grid!(df, coord.scratch_2d, coord, speed)
end

"""
    derivative!(df, f, coord, spectral)

Non-upwinding derivative.
"""
function derivative!(df, f, coord, spectral)
    # get the derivative at each grid point within each element and store in
    # coord.scratch_2d
    elementwise_derivative!(coord, f, spectral)
    # map the derivative from the elemental grid to the full grid;
    # at element boundaries, use the average of the derivatives from neighboring elements.
    derivative_elements_to_full_grid!(df, coord.scratch_2d, coord)
end

# Special versions for 'null' coordinates with only one point
function derivative!(df, f, coord, spectral::Union{null_spatial_dimension_info,
                                                   null_velocity_dimension_info})
    df .= 0.0
    return nothing
end

function second_derivative!(d2f, f, coord, spectral)
    # computes d^2f / d(coord)^2
    # For spectral element methods, calculate second derivative by applying first
    # derivative twice, with special treatment for element boundaries

    # First derivative
    elementwise_derivative!(coord, f, spectral)
    derivative_elements_to_full_grid!(coord.scratch3, coord.scratch_2d, coord)
    # MPI reconcile code here if used with z or r coords

    # Save elementwise first derivative result
    coord.scratch2_2d .= coord.scratch_2d

    # Second derivative for element interiors
    elementwise_derivative!(coord, coord.scratch3, spectral)
    derivative_elements_to_full_grid!(d2f, coord.scratch_2d, coord)
    # MPI reconcile code here if used with z or r coords

    if !isa(spectral, fourier_info)
        # The Fourier discretization does not have 'element boundaries' so does not need
        # this (attempted) 'correction'.
        # Add contribution to penalise discontinuous first derivatives at element
        # boundaries. For smooth functions this would do nothing so should not affect
        # convergence of the second derivative. Aims to stabilise numerical instability when
        # spike develops at an element boundary. The coefficient is an arbitrary choice, it
        # should probably be large enough for stability but as small as possible.
        #
        # Arbitrary numerical coefficient
        C = 1.0
        function penalise_discontinuous_first_derivative!(d2f, imin, imax, df)
            # Left element boundary
            d2f[imin] += C * df[1]

            # Right element boundary
            d2f[imax] -= C * df[end]

            return nothing
        end
        @views penalise_discontinuous_first_derivative!(d2f, 1, coord.imax[1],
                                                        coord.scratch2_2d[:,1])
        for ielement ∈ 2:coord.nelement_local
            @views penalise_discontinuous_first_derivative!(d2f, coord.imin[ielement]-1,
                                                            coord.imax[ielement],
                                                            coord.scratch2_2d[:,ielement])
        end

        if coord.periodic
            # Need to get first derivatives from opposite ends of grid
            if coord.nelement_local != coord.nelement_global
                error("Distributed memory MPI not yet supported here")
            end
            d2f[1] -= C * coord.scratch2_2d[end,end]
            # With the first derivative from the opposite end of the grid, d2f[end] here
            # should be equal to d2f[1] up to rounding errors...
            @debug_consistency_checks isapprox(d2f[end] + C * coord.scratch2_2d[1,1], d2f[1]; atol=1.0e-14)
            # ...but because arithmetic operations were in a different order, there may be
            # rounding errors, so set the two ends exactly equal to ensure consistency for the
            # rest of the code - we assume that duplicate versions of the 'same point' on
            # element boundaries (due to periodic bc or distributed-MPI block boundaries) are
            # exactly equal.
            d2f[end] = d2f[1]
        else
            # For stability don't contribute to evolution at boundaries, in case these
            # points are not set by a boundary condition.
            # Full grid may be across processes and bc only applied to extreme ends of the
            # domain.
            if coord.irank == 0
                d2f[1] = 0.0
            end
            if coord.irank == coord.nrank - 1
                d2f[end] = 0.0
            end
        end
    end
    return nothing
end

function laplacian_derivative!(d2f, f, coord, spectral)
    # computes (1/coord) d / coord ( coord d f / d(coord)) for vperp coordinate
    # For spectral element methods, calculate second derivative by applying first
    # derivative twice, with special treatment for element boundaries

    # First derivative
    elementwise_derivative!(coord, f, spectral)
    derivative_elements_to_full_grid!(coord.scratch3, coord.scratch_2d, coord)
    if coord.name == "vperp"
        # include the Jacobian
        @. coord.scratch3 *= coord.grid
    end
    # MPI reconcile code here if used with z or r coords

    # Save elementwise first derivative result
    coord.scratch2_2d .= coord.scratch_2d

    # Second derivative for element interiors
    elementwise_derivative!(coord, coord.scratch3, spectral)
    derivative_elements_to_full_grid!(d2f, coord.scratch_2d, coord)
    if coord.name == "vperp"
        # include the Jacobian
        @. d2f /= coord.grid
    end
    # MPI reconcile code here if used with z or r coords

    if !isa(spectral, fourier_info)
        # The Fourier discretization does not have 'element boundaries' so does not need
        # this (attempted) 'correction'.
        # Add contribution to penalise discontinuous first derivatives at element
        # boundaries. For smooth functions this would do nothing so should not affect
        # convergence of the second derivative. Aims to stabilise numerical instability when
        # spike develops at an element boundary. The coefficient is an arbitrary choice, it
        # should probably be large enough for stability but as small as possible.
        #
        # Arbitrary numerical coefficient
        C = 1.0
        function penalise_discontinuous_first_derivative!(d2f, imin, imax, df)
            # Left element boundary
            d2f[imin] += C * df[1]

            # Right element boundary
            d2f[imax] -= C * df[end]

            return nothing
        end
        @views penalise_discontinuous_first_derivative!(d2f, 1, coord.imax[1],
                                                        coord.scratch2_2d[:,1])
        for ielement ∈ 2:coord.nelement_local
            @views penalise_discontinuous_first_derivative!(d2f, coord.imin[ielement]-1,
                                                            coord.imax[ielement],
                                                            coord.scratch2_2d[:,ielement])
        end

        if coord.bc ∈ ("zero", "both_zero", "zero-no-regularity")
            # For stability don't contribute to evolution at boundaries, in case these
            # points are not set by a boundary condition.
            # Full grid may be across processes and bc only applied to extreme ends of the
            # domain.
            if coord.irank == coord.nrank - 1
                d2f[end] = 0.0
            end
        else
            error("Unsupported bc '$(coord.bc)'")
        end
    end
    return nothing
end

"""
A function that assigns the local weak-form matrices to a global array QQ_global for later
solving weak form of required 1D equation.

The 'option' variable is a flag for choosing the type of matrix to be constructed.
Currently the function is set up to assemble the elemental matrices without imposing
boundary conditions on the first and final rows of the matrix by default. This means that
the operators constructed from this function can only be used for differentiation, and not
solving 1D ODEs.  This assembly function assumes that the coordinate is not distributed.
To extend this function to support distributed-memory MPI, addition of off-memory matrix
elements to the exterior points would be required.

The typical use of this function is to assemble matrixes M and K in

 M * d2f = K * f

where M is the mass matrix and K is the stiffness matrix, and we wish to solve for d2f
given f. To solve 1D ODEs

K * f = b = M * d2f

for f given boundary data on f with periodic or dirichlet boundary conditions, set

`periodic_bc = true`, `b[end] = 0`

or

`dirichlet_bc = true`, `b[1] = f[1]` (except for cylindrical coordinates), `b[end] = f[end]`

in the function call, and create new matrices for this purpose in a
`weak_discretization_info` struct. Currently the Laplacian matrix is supported with
boundary conditions.
"""
function setup_global_weak_form_matrix!(QQ_global::Array{mk_float,2}, lobatto, radau,
                                        coord, option; dirichlet_bc=false,
                                        periodic_bc=false)
    QQ_j = allocate_float(coord.ngrid, coord.ngrid)

    ngrid = coord.ngrid
    imin = coord.imin
    imax = coord.imax
    @. QQ_global = 0.0

    # assembly below assumes no contributions from elements outside of local domain
    k = 0
    for j in 1:coord.nelement_local
        get_QQ_local!(QQ_j, j, lobatto, radau, coord, option)
        iminl = imin[j]-k
        imaxl = imax[j]
        @. QQ_global[iminl:imaxl,iminl:imaxl] += QQ_j[:,:]
        k = 1
    end

    if dirichlet_bc
        # Make matrix diagonal for first/last grid points so it does not change the values
        # there
        if !(coord.name == "vperp")
            # modify lower endpoint if not a radial/cylindrical coordinate
            if coord.irank == 0
                QQ_global[1,:] .= 0.0
                QQ_global[1,1] = 1.0
            end
        end
        if coord.irank == coord.nrank - 1
            QQ_global[end,:] .= 0.0
            QQ_global[end,end] = 1.0
        end
        # requires RHS vector b[1],b[end] = boundary values
    end
    if periodic_bc
        # Make periodic boundary condition by modifying elements of matrix for duplicate point
        # add assembly contribution to lower endpoint from upper endpoint
        j = coord.nelement_local
        get_QQ_local!(QQ_j, j, lobatto, radau, coord, option)
        iminl = imin[j] - mk_int(coord.nelement_local > 1)
        imaxl = imax[j]
        QQ_global[1,iminl:imaxl] .+= QQ_j[end,:]
        # Enforce continuity at the periodic boundary
        # All-zero row in RHS matrix sets last element of `b` vector in
        # `mass_matrix.x = b` to zero
        QQ_global[end,:] .= 0.0
        # enforce periodicity `x[1] = x[end]` using the last row of the
        # `A.x = b` system.
        QQ_global[end,1] = 1.0
        QQ_global[end,end] = -1.0
        if option == "L" # or any other derivative (ODE) matrix requiring two BC (periodicity + value at endpoint)
            QQ_global[1,:] .= 0.0
            QQ_global[1,1] = 1.0
        end
    end

    return nothing
end

"""
Construction function to provide the appropriate elemental matrix `Q` to the global matrix
assembly functions.
"""
function get_QQ_local!(QQ::AbstractArray{mk_float,2}, ielement, lobatto, radau, coord,
                       option)

    if option == "M"
        get_MM_local!(QQ, ielement, lobatto, radau, coord)
    elseif option == "R"
        get_MR_local!(QQ, ielement, lobatto, radau, coord)
    elseif option == "N"
        get_MN_local!(QQ, ielement, lobatto, radau, coord)
    elseif option == "P"
        get_PP_local!(QQ, ielement, lobatto, radau, coord)
    elseif option == "U"
        get_PU_local!(QQ, ielement, lobatto, radau, coord)
    elseif option == "S"
        get_SS_local!(QQ, ielement, lobatto, radau, coord)
    elseif option == "K"
        get_KK_local!(QQ, ielement, lobatto, radau, coord)
    elseif option == "K_with_BC_terms"
        get_KK_local!(QQ, ielement, lobatto, radau, coord, explicit_BC_terms=true)
    elseif option == "J"
        get_KJ_local!(QQ, ielement, lobatto, radau, coord)
    elseif option == "L"
        get_LL_local!(QQ, ielement, lobatto, radau, coord)
    elseif option == "L_with_BC_terms"
        get_LL_local!(QQ, ielement, lobatto, radau, coord, explicit_BC_terms=true)
    elseif option == "D"
        get_DD_local!(QQ, ielement, lobatto, radau, coord)
    end
    return nothing
end

"""
If called for `coord.name = vperp` elemental matrix `MM` on the \$i^{th}\$ element is
```math
 M_{jk} = \\int^{v_\\perp^U}_{v_\\perp^L}  \\varphi_j(v_\\perp)\\varphi_k(v_\\perp) v_\\perp d v_\\perp = \\int^1_{-1} (c_i + x s_i)l_j(x)l_k(x) s_i d x
```
with \$c_i\$ and \$s_i\$ the appropriate shift and scale factors, respectively.
Otherwise, if called for any other coordinate elemental matrix `MM` is
```math
 M_{jk} = \\int^{v_\\|^U}_{v_\\|^L}  \\varphi_j(v_\\|)\\varphi_k(v_\\|) d v_\\| = \\int^1_{-1} l_j(x)l_k(x) s_i d x.
```
"""
function get_MM_local!(QQ, ielement, lobatto, radau, coord)
    scale_factor = coord.element_scale[ielement]
    shift_factor = coord.element_shift[ielement]
    if coord.name == "vperp" # assume integrals of form int^infty_0 (.) vperp d vperp
        # extra scale and shift factors required because of vperp in integral
        if ielement > 1 || coord.irank > 0 # lobatto points
            @. QQ = (shift_factor*lobatto.M0 + scale_factor*lobatto.M1)*scale_factor
        else # radau points
            @. QQ = (shift_factor*radau.M0 + scale_factor*radau.M1)*scale_factor
        end
    else # assume integrals of form int^infty_-infty (.) d vpa
        @. QQ = lobatto.M0*scale_factor
    end
    return nothing
end

function get_SS_local!(QQ, ielement, lobatto, radau, coord)
    scale_factor = coord.element_scale[ielement]
    shift_factor = coord.element_shift[ielement]
    if coord.name == "vperp" # assume integrals of form int^infty_0 (.) vperp d vperp
        # extra scale and shift factors required because of vperp in integral
        if ielement > 1 || coord.irank > 0 # lobatto points
            @. QQ =  shift_factor*lobatto.S0 + scale_factor*lobatto.S1
            # boundary terms from integration by parts
            imin = coord.imin[ielement] - 1
            imax = coord.imax[ielement]
            QQ[1,1] -= coord.grid[imin]
            QQ[coord.ngrid,coord.ngrid] += coord.grid[imax]
        else # radau points
            @. QQ =  shift_factor*radau.S0 + scale_factor*radau.S1
            # boundary terms from integration by parts
            imax = coord.imax[ielement]
            QQ[coord.ngrid,coord.ngrid] += coord.grid[imax]
        end
    else # assume integrals of form int^infty_-infty (.) d vpa
        @. QQ = lobatto.S0
        # boundary terms from integration by parts
        QQ[1,1] -= 1.0
        QQ[coord.ngrid,coord.ngrid] += 1.0
    end
    return nothing
end

"""
If called for `coord.name = vperp` elemental matrix `KK` on the \$i^{th}\$ element is
```math
 K_{jk} = -\\int^{v_\\perp^U}_{v_\\perp^L} \\left(v_\\perp \\frac{\\partial\\varphi_j(v_\\perp)}{\\partial v_\\perp} + \\varphi_j(v_\\perp) \\right)
 \\frac{\\partial\\varphi_k(v_\\perp)}{\\partial v_\\perp} d v_\\perp
 = -\\int^1_{-1} ((c_i + x s_i)l_j^\\prime(x) + l_j(x))l_k^\\prime(x) d x /s_i
```
with \$c_i\$ and \$s_i\$ the appropriate shift and scale factors, respectively.
Otherwise, if called for any other coordinate elemental matrix `KK` is the same as `LL`
(see `get_LL_local!()).  If `explicit_BC_terms = true`, boundary terms arising from
integration by parts are included at the extreme boundary points.
"""
function get_KK_local!(QQ, ielement, lobatto, radau, coord; explicit_BC_terms=false)
    nelement = coord.nelement_local
    scale_factor = coord.element_scale[ielement]
    shift_factor = coord.element_shift[ielement]
    if coord.name == "vperp" # assume integrals of form int^infty_0 (.) vperp d vperp
        # extra scale and shift factors required because of vperp in integral
        # P0 factors make this a d^2 / dvperp^2 rather than (1/vperp) d ( vperp d (.) / d vperp)
        if ielement > 1 || coord.irank > 0 # lobatto points
            @. QQ =  (shift_factor/scale_factor)*lobatto.K0 + lobatto.K1 - lobatto.P0
            # boundary terms from integration by parts
            if explicit_BC_terms && ielement == 1 && coord.irank == 0
                imin = coord.imin[ielement] - 1
                @. QQ[1,:] -= coord.grid[imin]*lobatto.Dmat[1,:]/scale_factor
            end
            if explicit_BC_terms && ielement == nelement && coord.irank == coord.nrank - 1
                imax = coord.imax[ielement]
                @. QQ[coord.ngrid,:] += coord.grid[imax]*lobatto.Dmat[coord.ngrid,:]/scale_factor
            end
        else # radau points
            @. QQ =  (shift_factor/scale_factor)*radau.K0 + radau.K1 - radau.P0
            # boundary terms from integration by parts
            if explicit_BC_terms && ielement == nelement && coord.irank == coord.nrank - 1
                imax = coord.imax[ielement]
                @. QQ[coord.ngrid,:] += coord.grid[imax]*radau.Dmat[coord.ngrid,:]/scale_factor
            end
        end
    else # assume integrals of form int^infty_-infty (.) d vpa
        @. QQ = lobatto.K0/scale_factor
        if !coord.periodic
            # boundary terms from integration by parts
            if explicit_BC_terms && ielement == 1 && coord.irank == 0
                @. QQ[1,:] -= lobatto.Dmat[1,:]/scale_factor
            end
            if explicit_BC_terms && ielement == nelement && coord.irank == coord.nrank - 1
                @. QQ[coord.ngrid,:] += lobatto.Dmat[coord.ngrid,:]/scale_factor
            end
        end
    end
    return nothing
end

"""
If called for `coord.name = vperp` elemental matrix `KJ` on the \$i^{th}\$ element is
```math
 (KJ)_{jk} = -\\int^{v_\\perp^U}_{v_\\perp^L} \\frac{\\partial\\varphi_j(v_\\perp)}{\\partial v_\\perp}\\frac{\\partial\\varphi_k(v_\\perp)}{\\partial v_\\perp} v_\\perp^2 d v_\\perp
 = -\\int^1_{-1} (c_i + x s_i)^2l_j^\\prime(x)l_k^\\prime(x) d x /s_i
```
with \$c_i\$ and \$s_i\$ the appropriate shift and scale factors, respectively.
Otherwise, if called for any other coordinate elemental matrix `KJ` is the same as `LL`
(see `get_LL_local()!`).
"""
function get_KJ_local!(QQ, ielement, lobatto, radau, coord)
    scale_factor = coord.element_scale[ielement]
    shift_factor = coord.element_shift[ielement]
    if coord.name == "vperp" # assume integrals of form int^infty_0 (.) vperp d vperp
        # extra scale and shift factors required because of vperp^2 in integral
        if ielement > 1 || coord.irank > 0 # lobatto points
            @. QQ = (lobatto.K0*((shift_factor^2)/scale_factor) +
                     lobatto.K1*2.0*shift_factor +
                     lobatto.K2*scale_factor)
        else # radau points
            @. QQ =  (radau.K0*((shift_factor^2)/scale_factor) +
                      radau.K1*2.0*shift_factor +
                      radau.K2*scale_factor)
        end
    else # assume integrals of form int^infty_-infty (.) d vpa
        @. QQ = lobatto.K0/scale_factor
    end
    return nothing
end

"""
If called for `coord.name = vperp` elemental matrix `LL` on the \$i^{th}\$ element is
```math
 L_{jk} = -\\int^{v_\\perp^U}_{v_\\perp^L} \\frac{\\partial\\varphi_j(v_\\perp)}{\\partial v_\\perp}\\frac{\\partial\\varphi_k(v_\\perp)}{\\partial v_\\perp} v_\\perp d v_\\perp
 = -\\int^1_{-1} (c_i + x s_i)l_j^\\prime(x)l_k^\\prime(x) d x /s_i
```
with \$c_i\$ and \$s_i\$ the appropriate shift and scale factors, respectively.
Otherwise, if called for any other coordinate elemental matrix `LL` is
```math
 L_{jk} = -\\int^{v_\\|^U}_{v_\\|^L}  \\frac{\\partial\\varphi_j(v_\\|)}{\\partial v_\\|}\\frac{\\partial\\varphi_k(v_\\|)}{\\partial v_\\|} d v_\\| =
 -\\int^1_{-1} l_j^\\prime(x)l_k^\\prime(x) d x /s_i.
```
If `explicit_BC_terms = true`, boundary terms arising from integration by parts are
included at the extreme boundary points.
"""
function get_LL_local!(QQ, ielement, lobatto, radau, coord; explicit_BC_terms=false)
    nelement = coord.nelement_local
    scale_factor = coord.element_scale[ielement]
    shift_factor = coord.element_shift[ielement]
    if coord.name == "vperp" # assume integrals of form int^infty_0 (.) vperp d vperp
        # extra scale and shift factors required because of vperp in integral
        #  (1/vperp) d ( vperp d (.) / d vperp)
        if ielement > 1 || coord.irank > 0 # lobatto points
            @. QQ =  (shift_factor/scale_factor)*lobatto.K0 + lobatto.K1
            # boundary terms from integration by parts
            if explicit_BC_terms && ielement == 1 && coord.irank == 0
                imin = coord.imin[ielement] - 1
                @. QQ[1,:] -= coord.grid[imin]*lobatto.Dmat[1,:]/scale_factor
            end
            if explicit_BC_terms && ielement == nelement && coord.irank == coord.nrank - 1
                imax = coord.imax[ielement]
                @. QQ[coord.ngrid,:] += coord.grid[imax]*lobatto.Dmat[coord.ngrid,:]/scale_factor
            end
        else # radau points
            @. QQ =  (shift_factor/scale_factor)*radau.K0 + radau.K1
            # boundary terms from integration by parts
            if explicit_BC_terms && ielement == nelement && coord.irank == coord.nrank - 1
                imax = coord.imax[ielement]
                @. QQ[coord.ngrid,:] += coord.grid[imax]*radau.Dmat[coord.ngrid,:]/scale_factor
            end
        end
    else # d^2 (.) d vpa^2 -- assume integrals of form int^infty_-infty (.) d vpa
        @. QQ = lobatto.K0/scale_factor
        if !coord.periodic
            # boundary terms from integration by parts
            if explicit_BC_terms && ielement == 1 && coord.irank == 0
                @. QQ[1,:] -= lobatto.Dmat[1,:]/scale_factor
            end
            if explicit_BC_terms && ielement == nelement && coord.irank == coord.nrank - 1
                @. QQ[coord.ngrid,:] += lobatto.Dmat[coord.ngrid,:]/scale_factor
            end
        end
    end
    return nothing
end

# Strong-form differentiation matrix
function get_DD_local!(QQ, ielement, lobatto, radau, coord)
    scale_factor = coord.element_scale[ielement]
    if coord.name == "vperp" && ielement == 1 && coord.irank == 0
        @. QQ = radau.Dmat / scale_factor
    else
        @. QQ = lobatto.Dmat / scale_factor
    end
    return nothing
end

"""
If called for `coord.name = vperp` elemental matrix `MN` on the \$i^{th}\$ element is
```math
 (MN)_{jk} = \\int^{v_\\perp^U}_{v_\\perp^L}  \\varphi_j(v_\\perp)\\varphi_k(v_\\perp) d v_\\perp = \\int^1_{-1} l_j(x)l_k(x) s_i d x
```
with \$c_i\$ and \$s_i\$ the appropriate shift and scale factors, respectively.
Otherwise, if called for any other coordinate elemental matrix `MN` is the same as `MM`
(see `get_MM_local!()`).
"""
function get_MN_local!(QQ, ielement, lobatto, radau, coord)
    scale_factor = coord.element_scale[ielement]
    if coord.name == "vperp" # assume integrals of form int^infty_0 (.) vperp d vperp
        # extra scale and shift factors required because of vperp in integral
        if ielement > 1 || coord.irank > 0 # lobatto points
            @. QQ =  lobatto.M0*scale_factor
        else # radau points
            @. QQ =  radau.M0*scale_factor
        end
    else # assume integrals of form int^infty_-infty (.) d vpa
        @. QQ = lobatto.M0*scale_factor
    end
    return nothing
end

"""
If called for `coord.name = vperp` elemental matrix `MR` on the \$i^{th}\$ element is
```math
 (MR)_{jk} = \\int^{v_\\perp^U}_{v_\\perp^L}  \\varphi_j(v_\\perp)\\varphi_k(v_\\perp) v_\\perp^2 d v_\\perp = \\int^1_{-1} (c_i + s_i x)^2 l_j(x)l_k(x) s_i d x
```
with \$c_i\$ and \$s_i\$ the appropriate shift and scale factors, respectively.
Otherwise, if called for any other coordinate elemental matrix `MR` is the same as `MM`
(see `get_MM_local!()`).
"""
function get_MR_local!(QQ, ielement, lobatto, radau, coord)
    scale_factor = coord.element_scale[ielement]
    shift_factor = coord.element_shift[ielement]
    if coord.name == "vperp" # assume integrals of form int^infty_0 (.) vperp d vperp
        # extra scale and shift factors required because of vperp in integral
        if ielement > 1 || coord.irank > 0 # lobatto points
            @. QQ =  (lobatto.M0*shift_factor^2 +
                      lobatto.M1*2.0*shift_factor*scale_factor +
                      lobatto.M2*scale_factor^2)*scale_factor
        else # radau points
            @. QQ =  (radau.M0*shift_factor^2 +
                      radau.M1*2.0*shift_factor*scale_factor +
                      radau.M2*scale_factor^2)*scale_factor
        end
    else # assume integrals of form int^infty_-infty (.) d vpa
        @. QQ = lobatto.M0*scale_factor
    end
    return nothing
end

"""
If called for `coord.name = vperp` elemental matrix `PP` on the \$i^{th}\$ element is
```math
 P_{jk} = \\int^{v_\\perp^U}_{v_\\perp^L}  \\varphi_j(v_\\perp)\\frac{\\partial\\varphi_k(v_\\perp)}{\\partial v_\\perp} v_\\perp d v_\\perp
 = \\int^1_{-1} (c_i + x s_i)l_j(x)l_k^\\prime(x) d x
```
with \$c_i\$ and \$s_i\$ the appropriate shift and scale factors, respectively.
Otherwise, if called for any other coordinate elemental matrix `PP` is
```math
 P_{jk} = \\int^{v_\\|^U}_{v_\\|^L}  \\varphi_j(v_\\|)\\frac{\\partial\\varphi_k(v_\\|)}{\\partial v_\\|} d v_\\| = \\int^1_{-1} l_j(x)l_k^\\prime(x) d x.
```
"""
function get_PP_local!(QQ, ielement, lobatto, radau, coord)
    scale_factor = coord.element_scale[ielement]
    shift_factor = coord.element_shift[ielement]
    if coord.name == "vperp" # assume integrals of form int^infty_0 (.) vperp d vperp
        # extra scale and shift factors required because of vperp in integral
        if ielement > 1 || coord.irank > 0 # lobatto points
            @. QQ =  lobatto.P0*shift_factor + lobatto.P1*scale_factor
        else # radau points
            @. QQ =  radau.P0*shift_factor + radau.P1*scale_factor
        end
    else # assume integrals of form int^infty_-infty (.) d vpa
        @. QQ = lobatto.P0
    end
    return nothing
end

"""
If called for `coord.name = vperp` elemental matrix `PP` on the \$i^{th}\$ element is
```math
 (PU)_{jk} = \\int^{v_\\perp^U}_{v_\\perp^L}  \\varphi_j(v_\\perp)\\frac{\\partial\\varphi_k(v_\\perp)}{\\partial v_\\perp} v_\\perp^2 d v_\\perp
 = \\int^1_{-1} (c_i + x s_i)^2l_j(x)l_k^\\prime(x) d x
```
with \$c_i\$ and \$s_i\$ the appropriate shift and scale factors, respectively.
Otherwise, if called for any other coordinate elemental matrix `PU` is the same as `PP`
see `get_PP_local!()`.
"""
function get_PU_local!(QQ, ielement, lobatto, radau, coord)
    scale_factor = coord.element_scale[ielement]
    shift_factor = coord.element_shift[ielement]
    if coord.name == "vperp" # assume integrals of form int^infty_0 (.) vperp d vperp
        # extra scale and shift factors required because of vperp in integral
        if ielement > 1 || coord.irank > 0 # lobatto points
            @. QQ =  (lobatto.P0*shift_factor^2 +
                      lobatto.P1*2.0*shift_factor*scale_factor +
                      lobatto.P2*scale_factor^2)
        else # radau points
            @. QQ =  (radau.P0*shift_factor^2 +
                      radau.P1*2.0*shift_factor*scale_factor +
                      radau.P2*scale_factor^2)
        end
    else # assume integrals of form int^infty_-infty (.) d vpa
        @. QQ = lobatto.P0
    end
    return nothing
end

"""
Construction function for nonlinear diffusion matrices, only used in the assembly of the
collision operator
"""
function get_QQ_local!(QQ::AbstractArray{mk_float,3}, ielement, lobatto, radau, coord,
                       option)

    if option == "YY0" # mass-like matrix
        get_YY0_local!(QQ,ielement,lobatto,radau,coord)
    elseif option == "YY1" # first-derivative-like matrix
        get_YY1_local!(QQ,ielement,lobatto,radau,coord)
    elseif option == "YY2" # second-derivative-like matrix
        get_YY2_local!(QQ,ielement,lobatto,radau,coord)
    elseif option == "YY3" # first-derivative-like matrix
        get_YY3_local!(QQ,ielement,lobatto,radau,coord)
    end
    return nothing
end

"""
If called for `coord.name = vperp` elemental matrix `YY0` on the \$i^{th}\$ element is
```math
 (YY0)_{jkm} = \\int^{v_\\perp^U}_{v_\\perp^L}  \\varphi_j(v_\\perp)\\varphi_k(v_\\perp)\\varphi_m(v_\\perp) v_\\perp d v_\\perp
 = \\int^1_{-1} (c_i + x s_i)l_j(x)l_k(x)l_m(x) s_i d x
```
with \$c_i\$ and \$s_i\$ the appropriate shift and scale factors, respectively.
Otherwise, if called for any other coordinate elemental matrix `YY0` is
```math
 (YY0)_{jkm} = \\int^{v_\\|^U}_{v_\\|^L}  \\varphi_j(v_\\|)\\varphi_k(v_\\|)\\varphi_m(v_\\|) d v_\\|
 = \\int^1_{-1} l_j(x)l_k(x)l_m(x) s_i d x.
```
"""
function get_YY0_local!(QQ, ielement, lobatto, radau, coord)
    scale_factor = coord.element_scale[ielement]
    shift_factor = coord.element_shift[ielement]
    if coord.name == "vperp" # assume integrals of form int^infty_0 (.) vperp d vperp
        # extra scale and shift factors required because of vperp in integral
        if ielement > 1 || coord.irank > 0 # lobatto points
            @. QQ =  (shift_factor*lobatto.Y00 + scale_factor*lobatto.Y01)*scale_factor
        else # radau points
            @. QQ =  (shift_factor*radau.Y00 + scale_factor*radau.Y01)*scale_factor
        end
    else # assume integrals of form int^infty_-infty (.) d vpa
        @. QQ = lobatto.Y00*scale_factor
    end
    return nothing
end

"""
If called for `coord.name = vperp` elemental matrix `YY1` on the \$i^{th}\$ element is
```math
 (YY1)_{jkm} = \\int^{v_\\perp^U}_{v_\\perp^L}  \\varphi_j(v_\\perp)\\varphi_k(v_\\perp)\\frac{\\partial\\varphi_m(v_\\perp)}{\\partial v_\\perp} v_\\perp d v_\\perp
 = \\int^1_{-1} (c_i + x s_i)l_j(x)l_k(x)l_m^\\prime(x) d x
```
with \$c_i\$ and \$s_i\$ the appropriate shift and scale factors, respectively.
Otherwise, if called for any other coordinate elemental matrix `YY1` is
```math
 (YY1)_{jkm} = \\int^{v_\\|^U}_{v_\\|^L}  \\varphi_j(v_\\|)\\varphi_k(v_\\|)\\frac{\\partial\\varphi_m(v_\\|)}{\\partial v_\\|} d v_\\|
 = \\int^1_{-1} l_j(x)l_k(x)l_m^\\prime(x) d x.
```
"""
function get_YY1_local!(QQ, ielement, lobatto, radau, coord)
    scale_factor = coord.element_scale[ielement]
    shift_factor = coord.element_shift[ielement]
    if coord.name == "vperp" # assume integrals of form int^infty_0 (.) vperp d vperp
        # extra scale and shift factors required because of vperp in integral
        if ielement > 1 || coord.irank > 0 # lobatto points
            @. QQ =  shift_factor*lobatto.Y10 + scale_factor*lobatto.Y11
        else # radau points
            @. QQ =  shift_factor*radau.Y10 + scale_factor*radau.Y11
        end
    else # assume integrals of form int^infty_-infty (.) d vpa
        @. QQ = lobatto.Y10
    end
    return nothing
end

"""
If called for `coord.name = vperp` elemental matrix `YY2` on the \$i^{th}\$ element is
```math
 (YY2)_{jkm} = \\int^{v_\\perp^U}_{v_\\perp^L}  \\varphi_j(v_\\perp)\\frac{\\partial\\varphi_k(v_\\|)}{\\partial v_\\|}\\frac{\\partial\\varphi_m(v_\\perp)}{\\partial v_\\perp} v_\\perp d v_\\perp
 = \\int^1_{-1} (c_i + x s_i)l_j(x)l_k^\\prime(x)l_m^\\prime(x) d x/s_i
```
with \$c_i\$ and \$s_i\$ the appropriate shift and scale factors, respectively.
Otherwise, if called for any other coordinate elemental matrix `YY2` is
```math
 (YY2)_{jkm} = \\int^{v_\\|^U}_{v_\\|^L}  \\varphi_j(v_\\|)\\frac{\\partial\\varphi_k(v_\\|)}{\\partial v_\\|}\\frac{\\partial\\varphi_m(v_\\|)}{\\partial v_\\|} d v_\\|
 = \\int^1_{-1} l_j(x)l_k^\\prime(x)l_m^\\prime(x) d x /s_i.
```
"""
function get_YY2_local!(QQ, ielement, lobatto, radau, coord)
    scale_factor = coord.element_scale[ielement]
    shift_factor = coord.element_shift[ielement]
    if coord.name == "vperp" # assume integrals of form int^infty_0 (.) vperp d vperp
        # extra scale and shift factors required because of vperp in integral
        if ielement > 1 || coord.irank > 0 # lobatto points
            @. QQ =  (shift_factor/scale_factor)*lobatto.Y20 + lobatto.Y21
        else # radau points
            @. QQ =  (shift_factor/scale_factor)*radau.Y20 + radau.Y21
        end
    else # assume integrals of form int^infty_-infty (.) d vpa
        @. QQ = lobatto.Y20/scale_factor
    end
    return nothing
end

"""
If called for `coord.name = vperp` elemental matrix `YY3` on the \$i^{th}\$ element is
```math
 (YY3)_{jkm} = \\int^{v_\\perp^U}_{v_\\perp^L}  \\varphi_j(v_\\perp)\\frac{\\partial\\varphi_k(v_\\|)}{\\partial v_\\|}\\varphi_m(v_\\perp) v_\\perp d v_\\perp
 = \\int^1_{-1} (c_i + x s_i)l_j(x)l_k^\\prime(x)l_m(x) d x
```
with \$c_i\$ and \$s_i\$ the appropriate shift and scale factors, respectively.
Otherwise, if called for any other coordinate elemental matrix `YY3` is
```math
 (YY3)_{jkm} = \\int^{v_\\|^U}_{v_\\|^L}  \\varphi_j(v_\\|)\\frac{\\partial\\varphi_k(v_\\|)}{\\partial v_\\|}\\varphi_m(v_\\|) d v_\\|
 = \\int^1_{-1} l_j(x)l_k^\\prime(x)l_m(x) d x.
```
"""
function get_YY3_local!(QQ, ielement, lobatto, radau, coord)
    scale_factor = coord.element_scale[ielement]
    shift_factor = coord.element_shift[ielement]
    if coord.name == "vperp" # assume integrals of form int^infty_0 (.) vperp d vperp
        # extra scale and shift factors required because of vperp in integral
        if ielement > 1 || coord.irank > 0 # lobatto points
            @. QQ =  shift_factor*lobatto.Y30 + scale_factor*lobatto.Y31
        else # radau points
            @. QQ =  shift_factor*radau.Y30 + scale_factor*radau.Y31
        end
    else # assume integrals of form int^infty_-infty (.) d vpa
        @. QQ = lobatto.Y30
    end
    return nothing
end

"""
    mass_matrix_solve!(f, b, spectral::weak_discretization_info)

Solve
```math
M.f = b
```
for \$a\$, where \$M\$ is the mass matrix of a weak-form finite element method and \$b\$
is an input.
"""
function mass_matrix_solve!(f, b, spectral::weak_discretization_info)
    # invert mass matrix system
    ldiv!(f, spectral.mass_matrix_lu, b)
    return nothing
end

function second_derivative!(d2f, f, coord, spectral::weak_discretization_info)
    # obtain the RHS of numerical weak-form of the equation 
    # g = d^2 f / d coord^2, which is 
    # M * g = K * f, with M the mass matrix and K an appropriate stiffness matrix
    # by multiplying by basis functions and integrating by parts    
    mul!(coord.scratch3, spectral.K_matrix, f)

    # solve weak form matrix problem M * g = K * f to obtain g = d^2 f / d coord^2
    if coord.nrank > 1
        error("mass_matrix_solve!() does not support a "
              * "distributed coordinate")
    end
    # Do mass-matrix solve into a buffer array, to ensure that the output array is always
    # a contiguous array, not a view into another array that might have a stride bigger
    # than 1.
    mass_matrix_solve!(coord.scratch4, coord.scratch3, spectral)
    d2f .= coord.scratch4

    if coord.periodic
        # d2f[end] here should be equal to d2f[1] up to rounding errors...
        @debug_consistency_checks isapprox(d2f[end], d2f[1]; atol=1.0e-14)
        # ...but in the matrix operations arithmetic operations are not necessarily in
        # exactly the same order, there may be rounding errors, so set the two ends
        # exactly equal to ensure consistency for the rest of the code - we assume that
        # duplicate versions of the 'same point' on element boundaries (due to periodic bc
        # or distributed-MPI block boundaries) are exactly equal.
        d2f[end] = d2f[1]
    end
end

function laplacian_derivative!(d2f, f, coord, spectral::weak_discretization_info)
    # for coord.name 'vperp' obtain the RHS of numerical weak-form of the equation 
    # g = (1/coord) d/d coord ( coord  d f / d coord ), which is 
    # M * g = K * f, with M the mass matrix, and K an appropriate stiffness matrix,
    # by multiplying by basis functions and integrating by parts.
    # for all other coord.name, do exactly the same as second_derivative! above.
    mul!(coord.scratch3, spectral.L_matrix, f)

    if coord.periodic && coord.name == "vperp"
        error("laplacian_derivative!() cannot handle periodic boundaries for vperp")
    elseif coord.periodic
        if coord.nrank > 1
            error("laplacian_derivative!() cannot handle periodic boundaries for a "
                  * "distributed coordinate")
        end

        coord.scratch3[1] = 0.5 * (coord.scratch3[1] + coord.scratch3[end])
        coord.scratch3[end] = coord.scratch3[1]
    end

    # solve weak form matrix problem M * g = K * f to obtain g = d^2 f / d coord^2
    if coord.nrank > 1
        error("mass_matrix_solve!() does not support a "
              * "distributed coordinate")
    end
    mass_matrix_solve!(d2f, coord.scratch3, spectral)
end

"""
"""
function derivative_elements_to_full_grid!(df1d, df2d, coord, speed::AbstractArray{mk_float,1})
    # no changes need to be made for the derivative at points away from element boundaries
    elements_to_full_grid_interior_pts!(df1d, df2d, coord)
    # resolve the multi-valued nature of the derivative at element boundaries
    # by using the derivative from the upwind element
    reconcile_element_boundaries_upwind!(df1d, df2d, coord, speed)
    return nothing
end

"""
"""
function derivative_elements_to_full_grid!(df1d, df2d, coord)
    # no changes need to be made for the derivative at points away from element boundaries
    elements_to_full_grid_interior_pts!(df1d, df2d, coord)
    # resolve the multi-valued nature of the derivative at element boundaries
    # by using the derivative from the upwind element
    reconcile_element_boundaries_centered!(df1d, df2d, coord)
    return nothing
end

"""
maps the derivative at points away from element boundaries
from the grid/element representation to the full grid representation
"""
function elements_to_full_grid_interior_pts!(df1d, df2d, coord)
    # for efficiency, define ngm1 to be ngrid-1, as it will be used repeatedly
    ngm1 = coord.ngrid-1
    # treat the first element
    for i ∈ 2:ngm1
        df1d[i] = df2d[i,1]
    end
    # deal with any additional elements
    if coord.nelement_local > 1
        for ielem ∈ 2:coord.nelement_local
            for i ∈ 0:ngm1-2
                df1d[coord.imin[ielem]+i] = df2d[i+2,ielem]
            end
        end
    end
    return nothing
end

"""
if at the boundary point within the element, must carefully
choose which value of df to use; this is because
df is multi-valued at the overlapping point at the boundary
between neighboring elements.
here we choose to use the value of df from the upwind element.
"""
function reconcile_element_boundaries_upwind!(df1d, df2d, coord, speed::AbstractArray{mk_float,1})
    # note that the first ngrid points are classified as belonging to the first element
    # and the next ngrid-1 points belonging to second element, etc.

    # first deal with domain boundaries
    if coord.periodic && coord.nelement_global == coord.nelement_local
        # consider left domain boundary
        if speed[1] < 0.0
            # speed < 0 corresponds to negative advection speed, so
            # use derivative information from upwind element at larger coordinate value
            df1d[1] = df2d[1,1]
        elseif speed[1] > 0.0
            # speed > 0 corresponds to positive advection speed, so
            # use derivative information from upwind element at smaller coordinate value
            df1d[1] = df2d[coord.ngrid,coord.nelement_local]
        else
            # speed = 0, so no upwinding required;
            # use average value
            df1d[1] = 0.5*(df2d[1,1]+df2d[coord.ngrid,coord.nelement_local])
        end
        # consider right domain boundary
        if speed[coord.n] < 0.0
            # speed < 0 corresponds to negative advection speed, so
            # use derivative information from upwind element at larger coordinate value
            df1d[coord.n] = df2d[1,1]
        elseif speed[coord.n] > 0.0
            # speed > 0 corresponds to positive advection speed, so
            # use derivative information from upwind element at smaller coordinate value
            df1d[coord.n] = df2d[coord.ngrid,coord.nelement_local]
        else
            # speed = 0, so no upwinding required;
            # use average value
            df1d[coord.n] = 0.5*(df2d[1,1]+df2d[coord.ngrid,coord.nelement_local])
        end
    else
        df1d[1] = df2d[1,1]
        df1d[coord.n] = df2d[coord.ngrid,coord.nelement_local]
    end
    # next consider remaining elements, if any.
    # only need to consider interior element boundaries
    if coord.nelement_local > 1
        for ielem ∈ 2:coord.nelement_local
            im1 = ielem-1
            # consider left element boundary
            if speed[coord.imax[im1]] < 0.0
                # speed < 0 corresponds to negative advection speed, so
                # use derivative information from upwind element at larger coordinate value
                df1d[coord.imax[im1]] = df2d[1,ielem]
            elseif speed[coord.imax[im1]] > 0.0
                # speed > 0 corresponds to positive advection speed, so
                # use derivative information from upwind element at smaller coordinate value
                df1d[coord.imax[im1]] = df2d[coord.ngrid,im1]
            else
                # speed = 0, so no upwinding required;
                # use average value
                df1d[coord.imax[im1]] = 0.5*(df2d[1,ielem]+df2d[coord.ngrid,im1])
            end
        end
    end
    return nothing
end

"""
if at the boundary point within the element, must carefully
choose which value of df to use; this is because
df is multi-valued at the overlapping point at the boundary
between neighboring elements.
here we choose to use the averaged value across elements.
"""

function reconcile_element_boundaries_centered!(df1d, df2d, coord)
    # note that the first ngrid points are classified as belonging to the first element
    # and the next ngrid-1 points belonging to second element, etc.
	# first deal with domain boundaries
	if coord.periodic && coord.nelement_local == coord.nelement_global
		# consider left domain boundary
		df1d[1] = 0.5*(df2d[1,1]+df2d[coord.ngrid,coord.nelement_local])
		# consider right domain boundary
		df1d[coord.n] = df1d[1]
	else
	# put endpoints into 1D array to be reconciled
	# across processes at a higher scope -> larger message sizes possible
		df1d[1] = df2d[1,1]
		df1d[coord.n] = df2d[coord.ngrid,coord.nelement_local]
	end
	# next consider remaining elements, if any.
	# only need to consider interior element boundaries
	if coord.nelement_local > 1
		for ielem ∈ 2:coord.nelement_local
			im1 = ielem-1
			# consider left element boundary
			df1d[coord.imax[im1]] = 0.5*(df2d[1,ielem]+df2d[coord.ngrid,im1])
		end
	end
	return nothing
end

"""
extension of the above function to distributed memory MPI
function allows for arbitray array sizes ONLY IF the
if statements doing the final endpoint assignments are
updated to include each physical dimension required
in the main code
"""
function assign_endpoint!(df1d::AbstractArray{mk_float,Ndims},
                          receive_buffer::AbstractArray{mk_float,Mdims}, key::String,
                          coord; neutrals=false) where {Ndims,Mdims}
    if key == "lower"
            j = 1
    elseif key == "upper"
            j = coord.n
    else
            println("ERROR: invalid key in assign_endpoint!")
    end
    # test against coord name -- make sure to use exact string delimiters e.g. "x" not 'x'
    # test against Ndims (autodetermined) to choose which array slices to use in assigning endpoints
    if Ndims==1 && coord.name == "z"
        @begin_serial_region()
        @serial_region begin
            df1d[j] = receive_buffer[]
        end
    elseif Ndims==2 && coord.name == "z"
        @begin_r_region()
        @loop_r ir begin
            df1d[j,ir] = receive_buffer[ir]
        end
    elseif Ndims==3 && coord.name == "z"
        if neutrals
            @begin_sn_r_region()
            @loop_sn_r isn ir begin
                df1d[j,ir,isn] = receive_buffer[ir,isn]
            end
        else
            @begin_s_r_region()
            @loop_s_r is ir begin
                df1d[j,ir,is] = receive_buffer[ir,is]
            end
        end
    elseif Ndims==4 && coord.name == "z"
        @begin_r_vperp_vpa_region()
        @loop_r_vperp_vpa ir ivperp ivpa begin
            df1d[ivpa,ivperp,j,ir] = receive_buffer[ivpa,ivperp,ir]
        end
    elseif Ndims==5 && coord.name == "z"
        @begin_s_r_vperp_vpa_region()
        @loop_s_r_vperp_vpa is ir ivperp ivpa begin
            df1d[ivpa,ivperp,j,ir,is] = receive_buffer[ivpa,ivperp,ir,is]
        end
    elseif Ndims==6 && coord.name == "z"
        @begin_sn_r_vzeta_vr_vz_region()
        @loop_sn_r_vzeta_vr_vz isn ir ivzeta ivr ivz begin
            df1d[ivz,ivr,ivzeta,j,ir,isn] = receive_buffer[ivz,ivr,ivzeta,ir,isn]
        end
    elseif Ndims==1 && coord.name == "r"
        @begin_serial_region()
        @serial_region begin
            df1d[j] = receive_buffer[]
        end
    elseif Ndims==2 && coord.name == "r"
        @begin_z_region()
        @loop_z iz begin
            df1d[iz,j] = receive_buffer[iz]
        end
    elseif Ndims==3 && coord.name == "r"
        if neutrals
            @begin_sn_z_region()
            @loop_sn_z isn iz begin
                df1d[iz,j,isn] = receive_buffer[iz,isn]
            end
        else
            @begin_s_z_region()
            @loop_s_z is iz begin
                df1d[iz,j,is] = receive_buffer[iz,is]
            end
        end
    elseif Ndims==4 && coord.name == "r"
        @begin_z_vperp_vpa_region()
        @loop_z_vperp_vpa iz ivperp ivpa begin
            df1d[ivpa,ivperp,iz,j] = receive_buffer[ivpa,ivperp,iz]
        end
    elseif Ndims==5 && coord.name == "r"
        @begin_s_z_vperp_vpa_region()
        @loop_s_z_vperp_vpa is iz ivperp ivpa begin
            df1d[ivpa,ivperp,iz,j,is] = receive_buffer[ivpa,ivperp,iz,is]
        end
    elseif Ndims==6 && coord.name == "r"
        @begin_sn_z_vzeta_vr_vz_region()
        @loop_sn_z_vzeta_vr_vz isn iz ivzeta ivr ivz begin
            df1d[ivz,ivr,ivzeta,iz,j,isn] = receive_buffer[ivz,ivr,ivzeta,iz,isn]
        end
    else
        error("ERROR: failure to assign endpoints in reconcile_element_boundaries_MPI! (centered): coord.name: ",coord.name," Ndims: ",Ndims," key: ",key)
    end
end

function assign_endpoint_anyzv!(df1d::AbstractArray{mk_float,Ndims},
                                receive_buffer::AbstractArray{mk_float,Mdims},
                                key::String, coord) where {Ndims,Mdims}
    if key == "lower"
        j = 1
    elseif key == "upper"
        j = coord.n
    else
        println("ERROR: invalid key in assign_endpoint!")
    end
    # test against coord name -- make sure to use exact string delimiters e.g. "x" not 'x'
    # test against Ndims (autodetermined) to choose which array slices to use in assigning endpoints
    if Ndims==1 && coord.name == "z"
        @begin_anyzv_region()
        @anyzv_serial_region begin
            df1d[j] = receive_buffer[]
        end
    elseif Ndims==3 && coord.name == "z"
        @begin_anyzv_vperp_vpa_region()
        @loop_vperp_vpa ivperp ivpa begin
            df1d[ivpa,ivperp,j] = receive_buffer[ivpa,ivperp]
        end
    else
        error("ERROR: failure to assign endpoints in assign_endpoint_anyzv!: coord.name: ",coord.name," Ndims: ",Ndims," key: ",key)
    end
end

function apply_centred!(buffer::AbstractArray{mk_float,0},
                        endpoints::AbstractArray{mk_float,0}, coord_name::String;
                        neutrals=false)
    #buffer contains off-process endpoint
    #endpoint is local on-process endpoint
    #sgn = 1 for send irank -> irank + 1
    #sgn = -1 for send irank + 1 -> irank
    #loop over all indices in array
    @begin_serial_region()
    @serial_region begin
        buffer[] = 0.5*(buffer[] + endpoints[])
    end
end

function apply_centred!(buffer::AbstractArray{mk_float,1},
                        endpoints::AbstractArray{mk_float,1}, coord_name::String;
                        neutrals=false)
    #buffer contains off-process endpoint
    #endpoint is local on-process endpoint
    #sgn = 1 for send irank -> irank + 1
    #sgn = -1 for send irank + 1 -> irank
    #loop over all indices in array
    if coord_name == "r"
        @begin_z_region
        @loop_z iz begin
            buffer[iz] = 0.5*(buffer[iz] + endpoints[iz])
        end
    elseif coord_name == "z"
        @begin_r_region
        @loop_r ir begin
            buffer[ir] = 0.5*(buffer[ir] + endpoints[ir])
        end
    else
        error("Unsupported dimension '$coord_name'.")
    end
end

function apply_centred!(buffer::AbstractArray{mk_float,2},
                        endpoints::AbstractArray{mk_float,2}, coord_name::String;
                        neutrals)
    #buffer contains off-process endpoint
    #endpoint is local on-process endpoint
    #sgn = 1 for send irank -> irank + 1
    #sgn = -1 for send irank + 1 -> irank
    #loop over all indices in array
    if coord_name == "r"
        if neutrals
            @begin_sn_z_region
            @loop_sn_z isn iz begin
                buffer[iz,isn] = 0.5*(buffer[iz,isn] + endpoints[iz,isn])
            end
        else
            @begin_s_z_region
            @loop_s_z is iz begin
                buffer[iz,is] = 0.5*(buffer[iz,is] + endpoints[iz,is])
            end
        end
    elseif coord_name == "z"
        if neutrals
            @begin_sn_r_region
            @loop_sn_r isn ir begin
                buffer[ir,isn] = 0.5*(buffer[ir,isn] + endpoints[ir,isn])
            end
        else
            @begin_s_r_region
            @loop_s_r is ir begin
                buffer[ir,is] = 0.5*(buffer[ir,is] + endpoints[ir,is])
            end
        end
    else
        error("Unsupported dimension '$coord_name'.")
    end
end

function apply_centred!(buffer::AbstractArray{mk_float,3},
                        endpoints::AbstractArray{mk_float,3}, coord_name::String;
                        neutrals=false)
    #buffer contains off-process endpoint
    #endpoint is local on-process endpoint
    #sgn = 1 for send irank -> irank + 1
    #sgn = -1 for send irank + 1 -> irank
    #loop over all indices in array
    if coord_name == "r"
        @begin_z_vperp_vpa_region
        @loop_z_vperp_vpa iz ivperp ivpa begin
            buffer[ivpa,ivperp,iz] = 0.5*(buffer[ivpa,ivperp,iz] + endpoints[ivpa,ivperp,iz])
        end
    elseif coord_name == "z"
        @begin_r_vperp_vpa_region
        @loop_r_vperp_vpa ir ivperp ivpa begin
            buffer[ivpa,ivperp,ir] = 0.5*(buffer[ivpa,ivperp,ir] + endpoints[ivpa,ivperp,ir])
        end
    else
        error("Unsupported dimension '$coord_name'.")
    end
end

function apply_centred!(buffer::AbstractArray{mk_float,4},
                        endpoints::AbstractArray{mk_float,4}, coord_name::String;
                        neutrals=false)
    #buffer contains off-process endpoint
    #endpoint is local on-process endpoint
    #sgn = 1 for send irank -> irank + 1
    #sgn = -1 for send irank + 1 -> irank
    #loop over all indices in array
    if coord_name == "r"
        @begin_s_z_vperp_vpa_region
        @loop_s_z_vperp_vpa is iz ivperp ivpa begin
            buffer[ivpa,ivperp,iz,is] = 0.5*(buffer[ivpa,ivperp,iz,is] + endpoints[ivpa,ivperp,iz,is])
        end
    elseif coord_name == "z"
        @begin_s_r_vperp_vpa_region
        @loop_s_r_vperp_vpa is ir ivperp ivpa begin
            buffer[ivpa,ivperp,ir,is] = 0.5*(buffer[ivpa,ivperp,ir,is] + endpoints[ivpa,ivperp,ir,is])
        end
    else
        error("Unsupported dimension '$coord_name'.")
    end
end

function apply_centred!(buffer::AbstractArray{mk_float,5},
                        endpoints::AbstractArray{mk_float,5}, coord_name::String;
                        neutrals=false)
    #buffer contains off-process endpoint
    #endpoint is local on-process endpoint
    #sgn = 1 for send irank -> irank + 1
    #sgn = -1 for send irank + 1 -> irank
    #loop over all indices in array
    if coord_name == "r"
        @begin_sn_z_vzeta_vr_vz_region
        @loop_sn_z_vzeta_vr_vz isn iz ivzeta ivr ivz begin
            buffer[ivz,ivr,ivzeta,iz,isn] = 0.5*(buffer[ivz,ivr,ivzeta,iz,isn] + endpoints[ivz,ivr,ivzeta,iz,isn])
        end
    elseif coord_name == "z"
        @begin_sn_r_vzeta_vr_vz_region
        @loop_sn_r_vzeta_vr_vz isn ir ivzeta ivr ivz begin
            buffer[ivz,ivr,ivzeta,ir,isn] = 0.5*(buffer[ivz,ivr,ivzeta,ir,isn] + endpoints[ivz,ivr,ivzeta,ir,isn])
        end
    else
        error("Unsupported dimension '$coord_name'.")
    end
end

function apply_centred_anyzv!(buffer::AbstractArray{mk_float,0},
                              endpoints::AbstractArray{mk_float,0}, coord_name::String)
    #buffer contains off-process endpoint
    #endpoint is local on-process endpoint
    #sgn = 1 for send irank -> irank + 1
    #sgn = -1 for send irank + 1 -> irank
    #loop over all indices in array
    @begin_anyzv_region()
    @anyzv_serial_region begin
        buffer[] = 0.5*(buffer[] + endpoints[])
    end
end

function apply_centred_anyzv!(buffer::AbstractArray{mk_float,2},
                              endpoints::AbstractArray{mk_float,2}, coord_name::String)
    #buffer contains off-process endpoint
    #endpoint is local on-process endpoint
    #sgn = 1 for send irank -> irank + 1
    #sgn = -1 for send irank + 1 -> irank
    #loop over all indices in array
    if coord_name == "z"
        @begin_anyzv_vperp_vpa_region
        @loop_vperp_vpa ivperp ivpa begin
            buffer[ivpa,ivperp] = 0.5*(buffer[ivpa,ivperp] + endpoints[ivpa,ivperp])
        end
    else
        error("Unsupported dimension '$coord_name'.")
    end
end

@timeit_debug global_timer reconcile_element_boundaries_MPI!(
                  df1d::AbstractArray{mk_float,Ndims},
                  dfdx_lower_endpoints::AbstractArray{mk_float,Mdims},
                  dfdx_upper_endpoints::AbstractArray{mk_float,Mdims},
                  receive_buffer1::AbstractArray{mk_float,Mdims},
                  receive_buffer2::AbstractArray{mk_float,Mdims},
                  coord; neutrals=false) where {Ndims,Mdims} = begin
	
    nrank = coord.nrank
    irank = coord.irank

    # synchronize buffers
    # -- this all-to-all block communicate here requires that this function is NOT called from within a parallelised loop
    # -- or from a @serial_region or from an if statment isolating a single rank on a block
    @_block_synchronize()
    @serial_region begin

        # now deal with endpoints that are stored across ranks
        comm = coord.comm
        #send_buffer = coord.send_buffer
        #receive_buffer = coord.receive_buffer
        # sending pattern is cyclic. First we send data form irank -> irank + 1
        # to fix the lower endpoints, then we send data from irank -> irank - 1
        # to fix upper endpoints. Special exception for the periodic points.
        # receive_buffer[1] is for data received, send_buffer[1] is data to be sent

        # pass data from irank -> irank + 1, receive data from irank - 1
        rreq1 = MPI.Irecv!(receive_buffer1, comm; source=coord.prevrank, tag=1)
        sreq1 = MPI.Isend(dfdx_upper_endpoints, comm; dest=coord.nextrank, tag=1)

        # pass data from irank -> irank - 1, receive data from irank + 1
        rreq2 = MPI.Irecv!(receive_buffer2, comm; source=coord.nextrank, tag=2)
        sreq2 = MPI.Isend(dfdx_lower_endpoints, comm; dest=coord.prevrank, tag=2)
        stats = MPI.Waitall([rreq1, sreq1, rreq2, sreq2])
    end

    # synchronize buffers
    @_block_synchronize()

    # now update receive buffers, taking into account the reconciliation
    if irank == 0
        if coord.periodic
            #update the extreme lower endpoint with data from irank = nrank -1	
            apply_centred!(receive_buffer1, dfdx_lower_endpoints, coord.name; neutrals)
        else #directly use value from Cheb
            receive_buffer1 = dfdx_lower_endpoints
        end
    else # enforce continuity at lower endpoint
        apply_centred!(receive_buffer1, dfdx_lower_endpoints, coord.name; neutrals)
    end
    #now update the df1d array -- using a slice appropriate to the dimension reconciled
    assign_endpoint!(df1d, receive_buffer1, "lower", coord; neutrals)

    if irank == nrank-1
        if coord.periodic
            #update the extreme upper endpoint with data from irank = 0
            apply_centred!(receive_buffer2, dfdx_upper_endpoints, coord.name; neutrals)
        else #directly use value from Cheb
            receive_buffer2 = dfdx_upper_endpoints
        end
    else # enforce continuity at upper endpoint
        apply_centred!(receive_buffer2, dfdx_upper_endpoints, coord.name; neutrals)
    end
    #now update the df1d array -- using a slice appropriate to the dimension reconciled
    assign_endpoint!(df1d, receive_buffer2, "upper", coord; neutrals)
end

function apply_upwind!(buffer::AbstractArray{mk_float,0},
                       speed::AbstractArray{mk_float,0},
                       endpoints::AbstractArray{mk_float,0}, sgn::mk_int,
                       coord_name::String; neutrals=false)
    #buffer contains off-process endpoint
    #speed > 0 is positive advection speed
    #speed < 0 is negative advection speed
    #endpoint is local on-process endpoint
    #sgn = 1 for send irank -> irank + 1
    #sgn = -1 for send irank + 1 -> irank
    #loop over all indices in array
    @begin_serial_region
    @serial_region begin
        if sgn*speed[] < 0.0
            # replace buffer value with endpoint value
            buffer[] = endpoints[]
        elseif sgn*speed[] > 0.0
            #do nothing
        else #average values
            buffer[] = 0.5*(buffer[] + endpoints[])
        end
    end
end

@timeit_debug global_timer reconcile_element_boundaries_MPI_anyzv!(
                  df1d::AbstractArray{mk_float,Ndims},
                  dfdx_lower_endpoints::AbstractArray{mk_float,Mdims},
                  dfdx_upper_endpoints::AbstractArray{mk_float,Mdims},
                  receive_buffer1::AbstractArray{mk_float,Mdims},
                  receive_buffer2::AbstractArray{mk_float,Mdims},
                  coord) where {Ndims,Mdims} = begin

    # synchronize buffers
    # -- this all-to-all subblock communicate here requires that this function is only
    # -- called within an r-parallelised loop, not within z-, vperp- or vpa-loops or from
    # -- an @anyzv_serial_region or from an if statment isolating a single rank in a
    # -- subblock
    @_anyzv_subblock_synchronize()
    #if anyzv_subblock_rank[] == 0 # lead process on this shared-memory subblock
    @anyzv_serial_region begin

        # now deal with endpoints that are stored across ranks
        comm = coord.comm
        nrank = coord.nrank
        irank = coord.irank
        #send_buffer = coord.send_buffer
        #receive_buffer = coord.receive_buffer
        # sending pattern is cyclic. First we send data form irank -> irank + 1
        # to fix the lower endpoints, then we send data from irank -> irank - 1
        # to fix upper endpoints. Special exception for the periodic points.
        # receive_buffer[1] is for data received, send_buffer[1] is data to be sent

        # pass data from irank -> irank + 1, receive data from irank - 1
        rreq1 = MPI.Irecv!(receive_buffer1, comm; source=coord.prevrank, tag=1)
        sreq1 = MPI.Isend(dfdx_upper_endpoints, comm; dest=coord.nextrank, tag=1)

        # pass data from irank -> irank - 1, receive data from irank + 1
        rreq2 = MPI.Irecv!(receive_buffer2, comm; source=coord.nextrank, tag=2)
        sreq2 = MPI.Isend(dfdx_lower_endpoints, comm; dest=coord.prevrank, tag=2)
        stats = MPI.Waitall([rreq1, sreq1, rreq2, sreq2])

        # now update receive buffers, taking into account the reconciliation
        if irank == 0
            if coord.periodic
                #update the extreme lower endpoint with data from irank = nrank -1	
                apply_centred_anyzv!(receive_buffer1, dfdx_lower_endpoints, coord.name)
            else #directly use value from Cheb
                receive_buffer1 = dfdx_lower_endpoints
            end
        else # enforce continuity at lower endpoint
            apply_centred_anyzv!(receive_buffer1, dfdx_lower_endpoints, coord.name)
        end
        #now update the df1d array -- using a slice appropriate to the dimension reconciled
        assign_endpoint_anyzv!(df1d,receive_buffer1,"lower",coord)

        if irank == nrank-1
            if coord.periodic
                #update the extreme upper endpoint with data from irank = 0
                apply_centred_anyzv!(receive_buffer2, dfdx_upper_endpoints, coord.name)
            else #directly use value from Cheb
                receive_buffer2 = dfdx_upper_endpoints
            end
        else # enforce continuity at upper endpoint
            apply_centred_anyzv!(receive_buffer2, dfdx_upper_endpoints, coord.name)
        end
        #now update the df1d array -- using a slice appropriate to the dimension reconciled
        assign_endpoint_anyzv!(df1d,receive_buffer2,"upper",coord)

    end
    # synchronize buffers
    @_anyzv_subblock_synchronize()
end

function apply_upwind!(buffer::AbstractArray{mk_float,1},
                       speed::AbstractArray{mk_float,1},
                       endpoints::AbstractArray{mk_float,1}, sgn::mk_int,
                       coord_name::String; neutrals=false)
    #buffer contains off-process endpoint
    #speed > 0 is positive advection speed
    #speed < 0 is negative advection speed
    #endpoint is local on-process endpoint
    #sgn = 1 for send irank -> irank + 1
    #sgn = -1 for send irank + 1 -> irank
    #loop over all indices in array
    if coord_name == "r"
        @begin_z_region
        @loop_z iz begin
            if sgn*speed[iz] < 0.0
                # replace buffer value with endpoint value
                buffer[iz] = endpoints[iz]
            elseif sgn*speed[iz] > 0.0
                #do nothing
            else #average values
                buffer[iz] = 0.5*(buffer[iz] + endpoints[iz])
            end
        end
    elseif coord_name == "z"
        @begin_r_region
        @loop_r ir begin
            if sgn*speed[ir] < 0.0
                # replace buffer value with endpoint value
                buffer[ir] = endpoints[ir]
            elseif sgn*speed[ir] > 0.0
                #do nothing
            else #average values
                buffer[ir] = 0.5*(buffer[ir] + endpoints[ir])
            end
        end
    else
        error("Unsupported dimension '$coord_name'.")
    end
end

function apply_upwind!(buffer::AbstractArray{mk_float,2},
                       speed::AbstractArray{mk_float,2},
                       endpoints::AbstractArray{mk_float,2}, sgn::mk_int,
                       coord_name::String; neutrals)
    #buffer contains off-process endpoint
    #speed > 0 is positive advection speed
    #speed < 0 is negative advection speed
    #endpoint is local on-process endpoint
    #sgn = 1 for send irank -> irank + 1
    #sgn = -1 for send irank + 1 -> irank
    #loop over all indices in array
    if coord_name == "r"
        if neutrals
            @begin_sn_z_region
            @loop_sn_z isn iz begin
                if sgn*speed[iz,isn] < 0.0
                    # replace buffer value with endpoint value
                    buffer[iz,isn] = endpoints[iz,isn]
                elseif sgn*speed[iz,isn] > 0.0
                    #do nothing
                else #average values
                    buffer[iz,isn] = 0.5*(buffer[iz,isn] + endpoints[iz,isn])
                end
            end
        else
            @begin_s_z_region
            @loop_s_z is iz begin
                if sgn*speed[iz,is] < 0.0
                    # replace buffer value with endpoint value
                    buffer[iz,is] = endpoints[iz,is]
                elseif sgn*speed[iz,is] > 0.0
                    #do nothing
                else #average values
                    buffer[iz,is] = 0.5*(buffer[iz,is] + endpoints[iz,is])
                end
            end
        end
    elseif coord_name == "z"
        if neutrals
            @begin_sn_r_region
            @loop_sn_r isn ir begin
                if sgn*speed[ir,isn] < 0.0
                    # replace buffer value with endpoint value
                    buffer[ir,isn] = endpoints[ir,isn]
                elseif sgn*speed[ir,isn] > 0.0
                    #do nothing
                else #average values
                    buffer[ir,isn] = 0.5*(buffer[ir,isn] + endpoints[ir,isn])
                end
            end
        else
            @begin_s_r_region
            @loop_s_r is ir begin
                if sgn*speed[ir,is] < 0.0
                    # replace buffer value with endpoint value
                    buffer[ir,is] = endpoints[ir,is]
                elseif sgn*speed[ir,is] > 0.0
                    #do nothing
                else #average values
                    buffer[ir,is] = 0.5*(buffer[ir,is] + endpoints[ir,is])
                end
            end
        end
    else
        error("Unsupported dimension '$coord_name'.")
    end
end

function apply_upwind!(buffer::AbstractArray{mk_float,3},
                       speed::AbstractArray{mk_float,3},
                       endpoints::AbstractArray{mk_float,3}, sgn::mk_int,
                       coord_name::String; neutrals=false)
    #buffer contains off-process endpoint
    #speed > 0 is positive advection speed
    #speed < 0 is negative advection speed
    #endpoint is local on-process endpoint
    #sgn = 1 for send irank -> irank + 1
    #sgn = -1 for send irank + 1 -> irank
    #loop over all indices in array
    if coord_name == "r"
        @begin_z_vperp_vpa_region
        @loop_z_vperp_vpa iz ivperp ivpa begin
            if sgn*speed[ivpa,ivperp,iz] < 0.0
                # replace buffer value with endpoint value
                buffer[ivpa,ivperp,iz] = endpoints[ivpa,ivperp,iz]
            elseif sgn*speed[ivpa,ivperp,iz] > 0.0
                #do nothing
            else #average values
                buffer[ivpa,ivperp,iz] = 0.5*(buffer[ivpa,ivperp,iz] + endpoints[ivpa,ivperp,iz])
            end
        end
    elseif coord_name == "z"
        @begin_r_vperp_vpa_region
        @loop_r_vperp_vpa ir ivperp ivpa begin
            if sgn*speed[ivpa,ivperp,ir] < 0.0
                # replace buffer value with endpoint value
                buffer[ivpa,ivperp,ir] = endpoints[ivpa,ivperp,ir]
            elseif sgn*speed[ivpa,ivperp,ir] > 0.0
                #do nothing
            else #average values
                buffer[ivpa,ivperp,ir] = 0.5*(buffer[ivpa,ivperp,ir] + endpoints[ivpa,ivperp,ir])
            end
        end
    else
        error("Unsupported dimension '$coord_name'.")
    end
end

function apply_upwind!(buffer::AbstractArray{mk_float,4},
                       speed::AbstractArray{mk_float,4},
                       endpoints::AbstractArray{mk_float,4}, sgn::mk_int,
                       coord_name::String; neutrals=false)
    #buffer contains off-process endpoint
    #speed > 0 is positive advection speed
    #speed < 0 is negative advection speed
    #endpoint is local on-process endpoint
    #sgn = 1 for send irank -> irank + 1
    #sgn = -1 for send irank + 1 -> irank
    #loop over all indices in array
    if coord_name == "r"
        @begin_s_z_vperp_vpa_region
        @loop_s_z_vperp_vpa is iz ivperp ivpa begin
            if sgn*speed[ivpa,ivperp,iz,is] < 0.0
                # replace buffer value with endpoint value
                buffer[ivpa,ivperp,iz,is] = endpoints[ivpa,ivperp,iz,is]
            elseif sgn*speed[ivpa,ivperp,iz,is] > 0.0
                #do nothing
            else #average values
                buffer[ivpa,ivperp,iz,is] = 0.5*(buffer[ivpa,ivperp,iz,is] + endpoints[ivpa,ivperp,iz,is])
            end
        end
    elseif coord_name == "z"
        @begin_s_r_vperp_vpa_region
        @loop_s_r_vperp_vpa is ir ivperp ivpa begin
            if sgn*speed[ivpa,ivperp,ir,is] < 0.0
                # replace buffer value with endpoint value
                buffer[ivpa,ivperp,ir,is] = endpoints[ivpa,ivperp,ir,is]
            elseif sgn*speed[ivpa,ivperp,ir,is] > 0.0
                #do nothing
            else #average values
                buffer[ivpa,ivperp,ir,is] = 0.5*(buffer[ivpa,ivperp,ir,is] + endpoints[ivpa,ivperp,ir,is])
            end
        end
    else
        error("Unsupported dimension '$coord_name'.")
    end
end

function apply_upwind!(buffer::AbstractArray{mk_float,5},
                       speed::AbstractArray{mk_float,5},
                       endpoints::AbstractArray{mk_float,5}, sgn::mk_int,
                       coord_name::String; neutrals=false)
    #buffer contains off-process endpoint
    #speed > 0 is positive advection speed
    #speed < 0 is negative advection speed
    #endpoint is local on-process endpoint
    #sgn = 1 for send irank -> irank + 1
    #sgn = -1 for send irank + 1 -> irank
    #loop over all indices in array
    if coord_name == "r"
        @begin_sn_z_vzeta_vr_vz_region
        @loop_sn_z_vzeta_vr_vz isn iz ivzeta ivr ivz begin
            if sgn*speed[ivz,ivr,ivzeta,iz,isn] < 0.0
                # replace buffer value with endpoint value
                buffer[ivz,ivr,ivzeta,iz,isn] = endpoints[ivz,ivr,ivzeta,iz,isn]
            elseif sgn*speed[ivz,ivr,ivzeta,iz,isn] > 0.0
                #do nothing
            else #average values
                buffer[ivz,ivr,ivzeta,iz,isn] = 0.5*(buffer[ivz,ivr,ivzeta,iz,isn] + endpoints[ivz,ivr,ivzeta,iz,isn])
            end
        end
    elseif coord_name == "z"
        @begin_sn_r_vzeta_vr_vz_region
        @loop_sn_r_vzeta_vr_vz isn ir ivzeta ivr ivz begin
            if sgn*speed[ivz,ivr,ivzeta,ir,isn] < 0.0
                # replace buffer value with endpoint value
                buffer[ivz,ivr,ivzeta,ir,isn] = endpoints[ivz,ivr,ivzeta,ir,isn]
            elseif sgn*speed[ivz,ivr,ivzeta,ir,isn] > 0.0
                #do nothing
            else #average values
                buffer[ivz,ivr,ivzeta,ir,isn] = 0.5*(buffer[ivz,ivr,ivzeta,ir,isn] + endpoints[ivz,ivr,ivzeta,ir,isn])
            end
        end
    else
        error("Unsupported dimension '$coord_name'.")
    end
end

function apply_upwind_vpavperpz!(buffer::AbstractArray{mk_float,2},
                                 speed::AbstractArray{mk_float,2},
                                 endpoints::AbstractArray{mk_float,2}, sgn::mk_int)
    #buffer contains off-process endpoint
    #speed > 0 is positive advection speed
    #speed < 0 is negative advection speed
    #endpoint is local on-process endpoint
    #sgn = 1 for send irank -> irank + 1
    #sgn = -1 for send irank + 1 -> irank
    #loop over all indices in array
    @begin_anyzv_vperp_vpa_region()
    @loop_vperp_vpa ivperp ivpa begin
        if sgn*speed[ivpa,ivperp] < 0.0
            # replace buffer value with endpoint value
            buffer[ivpa,ivperp] = endpoints[ivpa,ivperp]
        elseif sgn*speed[ivpa,ivperp] > 0.0
            #do nothing
        else #average values
            buffer[ivpa,ivperp] = 0.5*(buffer[ivpa,ivperp] + endpoints[ivpa,ivperp])
        end
    end
    return nothing
end

@timeit_debug global_timer reconcile_element_boundaries_MPI!(
                  df1d::AbstractArray{mk_float,Ndims},
                  speed_lower_endpoints::AbstractArray{mk_float,Mdims},
                  speed_upper_endpoints::AbstractArray{mk_float,Mdims},
                  dfdx_lower_endpoints::AbstractArray{mk_float,Mdims},
                  dfdx_upper_endpoints::AbstractArray{mk_float,Mdims},
                  receive_buffer1::AbstractArray{mk_float,Mdims},
                  receive_buffer2::AbstractArray{mk_float,Mdims},
                  coord; neutrals=false) where {Ndims,Mdims} = begin
	
    nrank = coord.nrank
    irank = coord.irank

    # synchronize buffers
    # -- this all-to-all block communicate here requires that this function is NOT called from within a parallelised loop
    # -- or from a @serial_region or from an if statment isolating a single rank on a block
    @_block_synchronize()
    #if block_rank[] == 0 # lead process on this shared-memory block
    @serial_region begin
        # now deal with endpoints that are stored across ranks
        comm = coord.comm
        #send_buffer = coord.send_buffer
        #receive_buffer = coord.receive_buffer
        # sending pattern is cyclic. First we send data form irank -> irank + 1
        # to fix the lower endpoints, then we send data from irank -> irank - 1
        # to fix upper endpoints. Special exception for the periodic points.
        # receive_buffer[1] is for data received, send_buffer[1] is data to be sent

        # send highest end point on THIS rank
        # pass data from irank -> irank + 1, receive data from irank - 1
        rreq1 = MPI.Irecv!(receive_buffer1, comm; source=coord.prevrank, tag=1)
        sreq1 = MPI.Isend(dfdx_upper_endpoints, comm; dest=coord.nextrank, tag=1)

        # send lowest end point on THIS rank
        # pass data from irank -> irank - 1, receive data from irank + 1
        rreq2 = MPI.Irecv!(receive_buffer2, comm; source=coord.nextrank, tag=2)
        sreq2 = MPI.Isend(dfdx_lower_endpoints, comm; dest=coord.prevrank, tag=2)
        stats = MPI.Waitall([rreq1, sreq1, rreq2, sreq2])
    end

    # synchronize buffers
    @_block_synchronize()

    # now update receive buffers, taking into account the reconciliation
    if irank == 0
        if coord.periodic
            # depending on speed, update the extreme lower endpoint with data from irank = nrank -1	
            apply_upwind!(receive_buffer1, speed_lower_endpoints, dfdx_lower_endpoints, 1,
                          coord.name; neutrals)
        else # directly use value from Cheb at extreme lower point
            # Don't do an array copy here, just make the `receive_buffer1` variable refer
            # to a different array.
            receive_buffer1 = dfdx_lower_endpoints
        end
    else # depending on speed, update the lower endpoint with data from irank = nrank -1	
        apply_upwind!(receive_buffer1, speed_lower_endpoints, dfdx_lower_endpoints, 1,
                      coord.name; neutrals)
    end
    #now update the df1d array -- using a slice appropriate to the dimension reconciled
    assign_endpoint!(df1d, receive_buffer1, "lower", coord; neutrals)

    if irank == nrank-1
        if coord.periodic
            # depending on speed, update the extreme upper endpoint with data from irank = 0
            apply_upwind!(receive_buffer2, speed_upper_endpoints, dfdx_upper_endpoints,
                          -1, coord.name; neutrals)
        else #directly use value from Cheb
            # Don't do an array copy here, just make the `receive_buffer2` variable refer
            # to a different array.
            receive_buffer2 = dfdx_upper_endpoints
        end
    else # enforce continuity at upper endpoint
        apply_upwind!(receive_buffer2, speed_upper_endpoints, dfdx_upper_endpoints, -1,
                      coord.name; neutrals)
    end
    #now update the df1d array -- using a slice appropriate to the dimension reconciled
    assign_endpoint!(df1d, receive_buffer2, "upper", coord; neutrals)

    return nothing
end

# Special version for pdf_electron with no r-dimension, which has the same number of
# dimensions as an ion/neutral moment variable, but different dimensions.
@timeit_debug global_timer reconcile_element_boundaries_MPI_z_pdf_vpavperpz!(
                  df1d::AbstractArray{mk_float,3},
                  dfdx_lower_endpoints::AbstractArray{mk_float,2},
                  dfdx_upper_endpoints::AbstractArray{mk_float,2},
                  receive_buffer1::AbstractArray{mk_float,2},
                  receive_buffer2::AbstractArray{mk_float,2}, coord) = begin
	
    # synchronize buffers
    # -- this all-to-all subblock communicate here requires that this function is only
    # -- called within an r-parallelised loop, not within z-, vperp- or vpa-loops or from
    # -- an @anyzv_serial_region or from an if statment isolating a single rank in a
    # -- subblock
    @_anyzv_subblock_synchronize()
    #if anyzv_subblock_rank[] == 0 # lead process on this shared-memory subblock
    @anyzv_serial_region begin

        # now deal with endpoints that are stored across ranks
        comm = coord.comm
        nrank = coord.nrank
        irank = coord.irank
        #send_buffer = coord.send_buffer
        #receive_buffer = coord.receive_buffer
        # sending pattern is cyclic. First we send data form irank -> irank + 1
        # to fix the lower endpoints, then we send data from irank -> irank - 1
        # to fix upper endpoints. Special exception for the periodic points.
        # receive_buffer[1] is for data received, send_buffer[1] is data to be sent

        # pass data from irank -> irank + 1, receive data from irank - 1
        rreq1 = MPI.Irecv!(receive_buffer1, comm; source=coord.prevrank, tag=1)
        sreq1 = MPI.Isend(dfdx_upper_endpoints, comm; dest=coord.nextrank, tag=1)

        # pass data from irank -> irank - 1, receive data from irank + 1
        rreq2 = MPI.Irecv!(receive_buffer2, comm; source=coord.nextrank, tag=2)
        sreq2 = MPI.Isend(dfdx_lower_endpoints, comm; dest=coord.prevrank, tag=2)
        stats = MPI.Waitall([rreq1, sreq1, rreq2, sreq2])

        # now update receive buffers, taking into account the reconciliation
        if irank == 0
            if coord.periodic
                #update the extreme lower endpoint with data from irank = nrank -1	
                receive_buffer1 .= 0.5*(receive_buffer1 .+ dfdx_lower_endpoints)
            else #directly use value from Cheb
                receive_buffer1 .= dfdx_lower_endpoints
            end
        else # enforce continuity at lower endpoint
            receive_buffer1 .= 0.5*(receive_buffer1 .+ dfdx_lower_endpoints)
        end
        #now update the df1d array -- using a slice appropriate to the dimension reconciled
        @views df1d[:,:,1] .= receive_buffer1

        if irank == nrank-1
            if coord.periodic
                #update the extreme upper endpoint with data from irank = 0
                receive_buffer2 .= 0.5*(receive_buffer2 .+ dfdx_upper_endpoints)
            else #directly use value from Cheb
                receive_buffer2 .= dfdx_upper_endpoints
            end
        else # enforce continuity at upper endpoint
            receive_buffer2 .= 0.5*(receive_buffer2 .+ dfdx_upper_endpoints)
        end
        #now update the df1d array -- using a slice appropriate to the dimension reconciled
        @views df1d[:,:,end] .= receive_buffer2

    end
    # synchronize buffers
    @_anyzv_subblock_synchronize()
end

# Special version for pdf_electron with no r-dimension, which has the same number of
# dimensions as an ion/neutral moment variable, but different dimensions.
@timeit_debug global_timer reconcile_element_boundaries_MPI_z_pdf_vpavperpz!(
                  df1d::AbstractArray{mk_float,3},
                  speed_lower_endpoints::AbstractArray{mk_float,2},
                  speed_upper_endpoints::AbstractArray{mk_float,2},
                  dfdx_lower_endpoints::AbstractArray{mk_float,2},
                  dfdx_upper_endpoints::AbstractArray{mk_float,2},
                  receive_buffer1::AbstractArray{mk_float,2},
                  receive_buffer2::AbstractArray{mk_float,2}, coord) = begin
	
    # synchronize buffers
    # -- this all-to-all subblock communicate here requires that this function is only
    # -- called within an r-parallelised loop, not within z-, vperp- or vpa-loops or from
    # -- an @anyzv_serial_region or from an if statment isolating a single rank in a
    # -- subblock
    nrank = coord.nrank
    irank = coord.irank
    @begin_anyzv_region()
    @anyzv_serial_region begin
        # now deal with endpoints that are stored across ranks
        comm = coord.comm
        #send_buffer = coord.send_buffer
        #receive_buffer = coord.receive_buffer
        # sending pattern is cyclic. First we send data form irank -> irank + 1
        # to fix the lower endpoints, then we send data from irank -> irank - 1
        # to fix upper endpoints. Special exception for the periodic points.
        # receive_buffer[1] is for data received, send_buffer[1] is data to be sent

        # send highest end point on THIS rank
        # pass data from irank -> irank + 1, receive data from irank - 1
        rreq1 = MPI.Irecv!(receive_buffer1, comm; source=coord.prevrank, tag=1)
        sreq1 = MPI.Isend(dfdx_upper_endpoints, comm; dest=coord.nextrank, tag=1)

        # send lowest end point on THIS rank
        # pass data from irank -> irank - 1, receive data from irank + 1
        rreq2 = MPI.Irecv!(receive_buffer2, comm; source=coord.nextrank, tag=2)
        sreq2 = MPI.Isend(dfdx_lower_endpoints, comm; dest=coord.prevrank, tag=2)
        stats = MPI.Waitall([rreq1, sreq1, rreq2, sreq2])
    end

    # now update receive buffers, taking into account the reconciliation
    if irank == 0
        if coord.periodic
            # depending on speed, update the extreme lower endpoint with data from irank = nrank -1	
            apply_upwind_vpavperpz!(receive_buffer1, speed_lower_endpoints,
                                    dfdx_lower_endpoints, 1)
        else # directly use value from Cheb at extreme lower point
            # Don't do an array copy here, just make the `receive_buffer1` variable refer
            # to a different array.
            receive_buffer1 = dfdx_lower_endpoints
        end
    else # depending on speed, update the lower endpoint with data from irank = nrank -1	
        apply_upwind_vpavperpz!(receive_buffer1, speed_lower_endpoints,
                                dfdx_lower_endpoints, 1)
    end
    #now update the df1d array -- using a slice appropriate to the dimension reconciled
    @begin_anyzv_vperp_vpa_region()
    @loop_vperp_vpa ivperp ivpa begin
        df1d[ivpa,ivperp,1] = receive_buffer1[ivpa,ivperp]
    end

    if irank == nrank-1
        if coord.periodic
            # depending on speed, update the extreme upper endpoint with data from irank = 0
            apply_upwind_vpavperpz!(receive_buffer2, speed_upper_endpoints,
                                    dfdx_upper_endpoints, -1)
        else #directly use value from Cheb
            # Don't do an array copy here, just make the `receive_buffer2` variable refer
            # to a different array.
            receive_buffer2 = dfdx_upper_endpoints
        end
    else # enforce continuity at upper endpoint
        apply_upwind_vpavperpz!(receive_buffer2, speed_upper_endpoints,
                                dfdx_upper_endpoints, -1)
    end
    #now update the df1d array -- using a slice appropriate to the dimension reconciled
    @begin_anyzv_vperp_vpa_region()
    @loop_vperp_vpa ivperp ivpa begin
        df1d[ivpa,ivperp,end] = receive_buffer2[ivpa,ivperp]
    end

    # synchronize buffers
    @_anyzv_subblock_synchronize()
end

"""
Computes the integral of the integrand, using the input wgts
"""
function integral(integrand, wgts)
    # n is the number of grid points
    n = length(wgts)
    # initialize 'integral' to zero before sum
    integral = 0.0
    @debug_consistency_checks n == length(integrand) || throw(BoundsError(integrand))
    @debug_consistency_checks n == length(wgts) || throw(BoundsError(wgts))
    @inbounds for i ∈ 1:n
        integral += integrand[i]*wgts[i]
    end
    return integral
end

"""
Computes the integral of the integrand multiplied by v, using the input wgts
"""
function integral(integrand, v, wgts)
    # n is the number of grid points
    n = length(wgts)
    # initialize 'integral' to zero before sum
    integral = 0.0
    @debug_consistency_checks n == length(integrand) || throw(BoundsError(integrand))
    @debug_consistency_checks n == length(v) || throw(BoundsError(v))
    @debug_consistency_checks n == length(wgts) || throw(BoundsError(wgts))
    @inbounds for i ∈ 1:n
        integral += integrand[i] * v[i] * wgts[i]
    end
    return integral
end

"""
Computes the integral of the integrand multiplied by v^n, using the input wgts
"""
function integral(integrand, v, n, wgts)
    # n is the number of grid points
    n_v = length(wgts)
    # initialize 'integral' to zero before sum
    integral = 0.0
    @debug_consistency_checks n_v == length(integrand) || throw(BoundsError(integrand))
    @debug_consistency_checks n_v == length(v) || throw(BoundsError(v))
    @debug_consistency_checks n_v == length(wgts) || throw(BoundsError(wgts))
    @inbounds for i ∈ 1:n_v
        integral += integrand[i] * v[i] ^ n * wgts[i]
    end
    return integral
end

"""
Compute the 1D integral `∫dv prefactor(v)*integrand`

In this variant `v` should be a `coordinate` object.
"""
function integral(prefactor::Function, integrand, v)
    @debug_consistency_checks v.n == length(integrand) || throw(BoundsError(integrand))
    v_grid = v.grid
    wgts = v.wgts
    integral = 0.0
    @inbounds for i ∈ eachindex(v_grid)
        integral += prefactor(v[i]) * integrand[i] * wgts[i]
    end
    return integral
end

"""
Compute the 2D integral `∫d^2vperp.dvpa prefactor(vperp,vpa)*integrand`

In this variant `vperp` and `vpa` should be `coordinate` objects.

Note that vperp_wgts contains the extra factor of vperp required for the
Jacobian.
"""
function integral(prefactor::Function, integrand, vperp, vpa)
    @debug_consistency_checks (vpa.n, vperp.n) == size(integrand) || throw(BoundsError(integrand))
    vperp_grid = vperp.grid
    vperp_wgts = vperp.wgts
    vpa_grid = vpa.grid
    vpa_wgts = vpa.wgts
    integral = 0.0
    for ivperp ∈ eachindex(vperp_grid), ivpa ∈ eachindex(vpa_grid)
        integral += prefactor(vperp_grid[ivperp], vpa_grid[ivpa]) *
                    integrand[ivpa, ivperp] * vperp_wgts[ivperp] * vpa_wgts[ivpa]
    end
    return integral
end

"""
Compute the 3D integral `∫dvzeta.dvr.dvz prefactor(vzeta,vr,vz)*integrand`

In this variant `vzeta`, `vr`, and `vz` should be `coordinate` objects.
"""
function integral(prefactor::Function, integrand, vzeta, vr, vz)
    @debug_consistency_checks (vz.n, vr.n, vzeta.n) == size(integrand) || throw(BoundsError(integrand))
    vzeta_grid = vzeta.grid
    vzeta_wgts = vzeta.wgts
    vr_grid = vr.grid
    vr_wgts = vr.wgts
    vz_grid = vz.grid
    vz_wgts = vz.wgts
    integral = 0.0
    for ivzeta ∈ eachindex(vzeta_grid), ivr ∈ eachindex(vr_grid), ivz ∈ eachindex(vz_grid)
        integral += prefactor(vzeta_grid[ivzeta], vr_grid[ivr], vz_grid[ivz]) *
                    integrand[ivz, ivr, ivzeta] * vzeta_wgts[ivzeta] * vr_wgts[ivr] *
                    vz_wgts[ivz]
    end
    return integral
end

"""
2D velocity integration routines
"""

"""
Computes the integral of the 2D integrand, using the input wgts
"""
function integral(integrand, vx, px, wgtsx, vy, py, wgtsy)
    # nx is the number of grid points
    nx = length(wgtsx)
    ny = length(wgtsy)
    # initialize 'integral' to zero before sum
    integral = 0.0
    @debug_consistency_checks nx == size(integrand,1) || throw(BoundsError(integrand))
    @debug_consistency_checks ny == size(integrand,2) || throw(BoundsError(integrand))
    @debug_consistency_checks nx == length(vx) || throw(BoundsError(vx))
    @debug_consistency_checks ny == length(vy) || throw(BoundsError(vy))
#    @debug_consistency_checks ny == length(wgtsy) || throw(BoundsError(wtgsy))
#    @debug_consistency_checks nx == length(wgtsx) || throw(BoundsError(wtgsx))

    @inbounds for j ∈ 1:ny
        @inbounds for i ∈ 1:nx
            integral += integrand[i,j] * (vx[i] ^ px) * (vy[j] ^ py) * wgtsx[i] * wgtsy[j]
        end
    end
    return integral
end


"""
3D velocity integration routines
"""

"""
Computes the integral of the 3D integrand, using the input wgts
"""
function integral(integrand, vx, px, wgtsx, vy, py, wgtsy, vz, pz, wgtsz)
    # nx is the number of grid points
    nx = length(wgtsx)
    ny = length(wgtsy)
    nz = length(wgtsz)
    # initialize 'integral' to zero before sum
    integral = 0.0
    @debug_consistency_checks nx == size(integrand,1) || throw(BoundsError(integrand))
    @debug_consistency_checks ny == size(integrand,2) || throw(BoundsError(integrand))
    @debug_consistency_checks nz == size(integrand,3) || throw(BoundsError(integrand))
    @debug_consistency_checks nx == length(vx) || throw(BoundsError(vx))
    @debug_consistency_checks ny == length(vy) || throw(BoundsError(vy))
    @debug_consistency_checks nz == length(vz) || throw(BoundsError(vz))

    @inbounds for k ∈ 1:nz
        @inbounds for j ∈ 1:ny
            @inbounds for i ∈ 1:nx
                integral += integrand[i,j,k] * (vx[i] ^ px) * (vy[j] ^ py) * (vz[k] ^ pz) * wgtsx[i] * wgtsy[j] * wgtsz[k]
            end
        end
    end
    return integral
end


# Indefinite integral routines

"""
A function that takes the indefinite integral in each element of `coord.grid`,
leaving the result (element-wise) in `coord.scratch_2d`.
"""
function elementwise_indefinite_integration!(coord, ff, spectral::discretization_info)
    # the primitive of f
    pf = coord.scratch_2d
    # define local variable nelement for convenience
    nelement = coord.nelement_local
    # check array bounds
    @debug_consistency_checks nelement == size(pf,2) && coord.ngrid == size(pf,1) || throw(BoundsError(pf))

    # variable k will be used to avoid double counting of overlapping point
    k = 0
    j = 1 # the first element
    imin = coord.imin[j]-k
    # imax is the maximum index on the full grid for this (jth) element
    imax = coord.imax[j]
    if coord.radau_first_element && coord.irank == 0 # differentiate this element with the Radau scheme
        @views mul!(pf[:,j],spectral.radau.indefinite_integration_matrix[:,:],ff[imin:imax])
    else #differentiate using the Lobatto scheme
        @views mul!(pf[:,j],spectral.lobatto.indefinite_integration_matrix[:,:],ff[imin:imax])
    end
    # transform back to the physical coordinate scale
    for i in 1:coord.ngrid
        pf[i,j] *= coord.element_scale[j]
    end
    # calculate the derivative on each element
    @inbounds for j ∈ 2:nelement
        k = 1
        imin = coord.imin[j]-k
        # imax is the maximum index on the full grid for this (jth) element
        imax = coord.imax[j]
        @views mul!(pf[:,j],spectral.lobatto.indefinite_integration_matrix[:,:],ff[imin:imax])
        # transform back to the physical coordinate scale
        for i in 1:coord.ngrid
            pf[i,j] *= coord.element_scale[j]
        end
    end
    return nothing
end

end
