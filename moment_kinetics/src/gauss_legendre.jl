"""
module for Gauss-Legendre-Lobatto and Gauss-Legendre-Radau spectral element grids
"""
module gauss_legendre

export gausslobattolegendre_differentiation_matrix!
export gaussradaulegendre_differentiation_matrix!
export GaussLegendreLobatto_mass_matrix!
export GaussLegendre_mass_matrix_1!
export GaussLegendreLobatto_inverse_mass_matrix!
export GaussLegendreLobatto_K_matrix!
export scaled_gauss_legendre_lobatto_grid
export scaled_gauss_legendre_radau_grid
export gausslegendre_derivative!
export gausslegendre_apply_Kmat!
export gausslegendre_apply_Lmat!
export setup_gausslegendre_pseudospectral
export GaussLegendre_weak_product_matrix!
export ielement_global_func
export get_QQ_local!
export integration_matrix!

using FastGaussQuadrature
using LegendrePolynomials: Pl, dnPl
using LinearAlgebra: mul!, lu, ldiv!
using SparseArrays: sparse, AbstractSparseArray
using SparseMatricesCSR
using ..type_definitions: mk_float, mk_int
using ..array_allocation: allocate_float
import ..calculus: elementwise_derivative!, mass_matrix_solve!,
                   elementwise_indefinite_integration!,
                   setup_global_weak_form_matrix!
using ..debugging
import ..interpolation: single_element_interpolate!,
                        fill_single_element_interpolation_matrix!
using ..lagrange_polynomials: lagrange_poly_optimised, lagrange_poly_derivative_optimised, lagrange_poly
using ..moment_kinetics_structs: weak_discretization_info


#structs for passing around matrices for taking
#the derivatives on Gauss-Legendre points in 1D
"""
A struct for passing around elemental matrices
on Gauss-Legendre points in 1D
"""
struct gausslegendre_base_info
    # flag for whether we are Radau or Lobatto
    is_lobatto::Bool
    # elementwise differentiation matrix (ngrid*ngrid)
    Dmat::Array{mk_float,2}
    # elementwise integration matrix (ngrid*ngrid)
    indefinite_integration_matrix::Array{mk_float,2}
    # local mass matrix type 0
    M0::Array{mk_float,2}
    # local mass matrix type 1
    M1::Array{mk_float,2}
    # local mass matrix type 2
    M2::Array{mk_float,2}
    # local S (weak derivative) matrix type 0
    S0::Array{mk_float,2}
    # local S (weak derivative) matrix type 1
    S1::Array{mk_float,2}
    # local K (weak second derivative) matrix type 0
    K0::Array{mk_float,2}
    # local K (weak second derivative) matrix type 1
    K1::Array{mk_float,2}
    # local K (weak second derivative) matrix type 2
    K2::Array{mk_float,2}
    # local P (weak derivative no integration by parts) matrix type 0
    P0::Array{mk_float,2}
    # local P (weak derivative no integration by parts) matrix type 1
    P1::Array{mk_float,2}
    # local P (weak derivative no integration by parts) matrix type 2
    P2::Array{mk_float,2}
    # boundary condition differentiation matrix (for vperp grid using radau points)
    D0::Array{mk_float,1}
    # local nonlinear diffusion matrix Y00
    Y00::Array{mk_float,3}
    # local nonlinear diffusion matrix Y01
    Y01::Array{mk_float,3}
    # local nonlinear diffusion matrix Y10
    Y10::Array{mk_float,3}
    # local nonlinear diffusion matrix Y11
    Y11::Array{mk_float,3}
    # local nonlinear diffusion matrix Y20
    Y20::Array{mk_float,3}
    # local nonlinear diffusion matrix Y21
    Y21::Array{mk_float,3}
    # local nonlinear diffusion matrix Y30
    Y30::Array{mk_float,3}
    # local nonlinear diffusion matrix Y31
    Y31::Array{mk_float,3}
end

"""
A struct for Gauss-Legendre arrays needed for global operations in 1D,
contains the struct of elemental matrices for Lobatto and Radau points,
as well as some assembled 1D global matrices.
"""
struct gausslegendre_info{TSparse, TLU} <: weak_discretization_info
    lobatto::gausslegendre_base_info
    radau::gausslegendre_base_info
    # global (1D) mass matrix
    mass_matrix::TSparse
    # global (1D) weak second derivative matrix
    K_matrix::TSparse
    # global (1D) weak Laplacian derivative matrix
    L_matrix::TSparse
    # global (1D) weak second derivative matrix, with inverse mass matrix included (so
    # matrix is dense)
    dense_second_deriv_matrix::Array{mk_float,2}
    # mass matrix global (1D) LU object
    mass_matrix_lu::TLU
end

"""
Function to create `gausslegendre_info` struct.
"""
function setup_gausslegendre_pseudospectral(coord; collision_operator_dim=true)
    lobatto = setup_gausslegendre_pseudospectral_lobatto(coord,collision_operator_dim=collision_operator_dim)
    radau = setup_gausslegendre_pseudospectral_radau(coord,collision_operator_dim=collision_operator_dim)

    mass_matrix = allocate_float(coord, coord)
    K_matrix = allocate_float(coord, coord)
    L_matrix = allocate_float(coord, coord)

    dirichlet_bc = (coord.bc in ["zero", "both_zero", "constant"]) # and further options in future
    periodic_bc = coord.periodic
    setup_global_weak_form_matrix!(mass_matrix, lobatto, radau, coord, "M"; periodic_bc=periodic_bc)
    setup_global_weak_form_matrix!(K_matrix, lobatto, radau, coord, "K_with_BC_terms"; periodic_bc=periodic_bc)
    setup_global_weak_form_matrix!(L_matrix, lobatto, radau, coord, "L_with_BC_terms")
    dense_second_deriv_matrix = mass_matrix \ K_matrix

    mass_matrix = sparse(mass_matrix)
    mass_matrix_lu = lu(sparse(mass_matrix))

    return gausslegendre_info(lobatto, radau, mass_matrix, sparse(K_matrix),
                              sparse(L_matrix), dense_second_deriv_matrix, mass_matrix_lu)
end

"""
Function that fills and `n` x `n` array with the values of the identity matrix `I`.
"""
function identity_matrix!(I,n)
    @. I[:,:] = 0.0
    for i in 1:n
       I[i,i] = 1.0
    end
    return nothing    
end

"""
Function that creates the `gausslegendre_base_info` struct for Lobatto points.
If `collision_operator_dim = true`, assign the elemental matrices used to
implement the Fokker-Planck collision operator.
"""
function setup_gausslegendre_pseudospectral_lobatto(coord; collision_operator_dim=true)
    x, w = gausslobatto(coord.ngrid)
    x = mk_float.(x)
    w = mk_float.(w)
    Dmat = allocate_float(coord.ngrid, coord.ngrid)
    gausslobattolegendre_differentiation_matrix!(Dmat,x,coord.ngrid)
    indefinite_integration_matrix = allocate_float(coord.ngrid, coord.ngrid)
    integration_matrix!(indefinite_integration_matrix,x,coord.ngrid)
    
    M0 = allocate_float(coord.ngrid, coord.ngrid)
    GaussLegendre_weak_product_matrix!(M0,coord.ngrid,x,w,"M0")
    M1 = allocate_float(coord.ngrid, coord.ngrid)
    GaussLegendre_weak_product_matrix!(M1,coord.ngrid,x,w,"M1")
    K0 = allocate_float(coord.ngrid, coord.ngrid)
    GaussLegendre_weak_product_matrix!(K0,coord.ngrid,x,w,"K0")
    K1 = allocate_float(coord.ngrid, coord.ngrid)
    GaussLegendre_weak_product_matrix!(K1,coord.ngrid,x,w,"K1")
    P0 = allocate_float(coord.ngrid, coord.ngrid)
    GaussLegendre_weak_product_matrix!(P0,coord.ngrid,x,w,"P0")
    D0 = allocate_float(coord.ngrid)
    #@. D0 = Dmat[1,:] # values at lower extreme of element
    GaussLegendre_derivative_vector!(D0,-1.0,coord.ngrid,x,w)
    
    if collision_operator_dim
        M2 = allocate_float(coord.ngrid, coord.ngrid)
        GaussLegendre_weak_product_matrix!(M2,coord.ngrid,x,w,"M2")
        S0 = allocate_float(coord.ngrid, coord.ngrid)
        GaussLegendre_weak_product_matrix!(S0,coord.ngrid,x,w,"S0")
        S1 = allocate_float(coord.ngrid, coord.ngrid)
        GaussLegendre_weak_product_matrix!(S1,coord.ngrid,x,w,"S1")
        K2 = allocate_float(coord.ngrid, coord.ngrid)
        GaussLegendre_weak_product_matrix!(K2,coord.ngrid,x,w,"K2")
        P1 = allocate_float(coord.ngrid, coord.ngrid)
        GaussLegendre_weak_product_matrix!(P1,coord.ngrid,x,w,"P1")
        P2 = allocate_float(coord.ngrid, coord.ngrid)
        GaussLegendre_weak_product_matrix!(P2,coord.ngrid,x,w,"P2")

        Y00 = allocate_float(coord.ngrid, coord.ngrid, coord.ngrid)
        Y01 = allocate_float(coord.ngrid, coord.ngrid, coord.ngrid)
        Y10 = allocate_float(coord.ngrid, coord.ngrid, coord.ngrid)
        Y11 = allocate_float(coord.ngrid, coord.ngrid, coord.ngrid)
        Y20 = allocate_float(coord.ngrid, coord.ngrid, coord.ngrid)
        Y21 = allocate_float(coord.ngrid, coord.ngrid, coord.ngrid)
        Y30 = allocate_float(coord.ngrid, coord.ngrid, coord.ngrid)
        Y31 = allocate_float(coord.ngrid, coord.ngrid, coord.ngrid)
        GaussLegendre_weak_product_matrix!(Y00,coord.ngrid,x,w,"Y00")
        GaussLegendre_weak_product_matrix!(Y01,coord.ngrid,x,w,"Y01")
        GaussLegendre_weak_product_matrix!(Y10,coord.ngrid,x,w,"Y10")
        GaussLegendre_weak_product_matrix!(Y11,coord.ngrid,x,w,"Y11")
        GaussLegendre_weak_product_matrix!(Y20,coord.ngrid,x,w,"Y20")
        GaussLegendre_weak_product_matrix!(Y21,coord.ngrid,x,w,"Y21")
        GaussLegendre_weak_product_matrix!(Y30,coord.ngrid,x,w,"Y30")
        GaussLegendre_weak_product_matrix!(Y31,coord.ngrid,x,w,"Y31")
    else
        M2 = allocate_float(0, 0)
        S0 = allocate_float(0, 0)
        S1 = allocate_float(0, 0)
        K2 = allocate_float(0, 0)
        P1 = allocate_float(0, 0)
        P2 = allocate_float(0, 0)

        Y00 = allocate_float(0, 0, 0)
        Y01 = allocate_float(0, 0, 0)
        Y10 = allocate_float(0, 0, 0)
        Y11 = allocate_float(0, 0, 0)
        Y20 = allocate_float(0, 0, 0)
        Y21 = allocate_float(0, 0, 0)
        Y30 = allocate_float(0, 0, 0)
        Y31 = allocate_float(0, 0, 0)
    end
    return gausslegendre_base_info(true,Dmat,indefinite_integration_matrix, M0, M1, M2, S0, S1,
                                   K0, K1, K2, P0, P1, P2, D0, Y00, Y01, Y10, Y11, Y20,
                                   Y21, Y30, Y31)
end

"""
Function that creates the `gausslegendre_base_info` struct for Lobatto points.
If `collision_operator_dim = true`, assign the elemental matrices used to
implement the Fokker-Planck collision operator.
"""
function setup_gausslegendre_pseudospectral_radau(coord; collision_operator_dim=true)
    # Gauss-Radau points on [-1,1)
    x, w = gaussradau(coord.ngrid)
    x = mk_float.(x)
    w = mk_float.(w)
    # Gauss-Radau points on (-1,1] 
    xreverse, wreverse = -reverse(x), reverse(w)
    # elemental differentiation matrix
    Dmat = allocate_float(coord.ngrid, coord.ngrid)
    gaussradaulegendre_differentiation_matrix!(Dmat,x,coord.ngrid)
    indefinite_integration_matrix = allocate_float(coord.ngrid, coord.ngrid)
    integration_matrix!(indefinite_integration_matrix,xreverse,coord.ngrid)

    M0 = allocate_float(coord.ngrid, coord.ngrid)
    GaussLegendre_weak_product_matrix!(M0,coord.ngrid,xreverse,wreverse,"M0",radau=true)
    M1 = allocate_float(coord.ngrid, coord.ngrid)
    GaussLegendre_weak_product_matrix!(M1,coord.ngrid,xreverse,wreverse,"M1",radau=true)
    K0 = allocate_float(coord.ngrid, coord.ngrid)
    GaussLegendre_weak_product_matrix!(K0,coord.ngrid,xreverse,wreverse,"K0",radau=true)
    K1 = allocate_float(coord.ngrid, coord.ngrid)
    GaussLegendre_weak_product_matrix!(K1,coord.ngrid,xreverse,wreverse,"K1",radau=true)
    P0 = allocate_float(coord.ngrid, coord.ngrid)
    GaussLegendre_weak_product_matrix!(P0,coord.ngrid,xreverse,wreverse,"P0",radau=true)
    D0 = allocate_float(coord.ngrid)
    GaussLegendre_derivative_vector!(D0,-1.0,coord.ngrid,xreverse,wreverse,radau=true)
    if collision_operator_dim
        M2 = allocate_float(coord.ngrid, coord.ngrid)
        GaussLegendre_weak_product_matrix!(M2,coord.ngrid,xreverse,wreverse,"M2",radau=true)
        S0 = allocate_float(coord.ngrid, coord.ngrid)
        GaussLegendre_weak_product_matrix!(S0,coord.ngrid,xreverse,wreverse,"S0",radau=true)
        S1 = allocate_float(coord.ngrid, coord.ngrid)
        GaussLegendre_weak_product_matrix!(S1,coord.ngrid,xreverse,wreverse,"S1",radau=true)
        K2 = allocate_float(coord.ngrid, coord.ngrid)
        GaussLegendre_weak_product_matrix!(K2,coord.ngrid,xreverse,wreverse,"K2",radau=true)
        P1 = allocate_float(coord.ngrid, coord.ngrid)
        GaussLegendre_weak_product_matrix!(P1,coord.ngrid,xreverse,wreverse,"P1",radau=true)
        P2 = allocate_float(coord.ngrid, coord.ngrid)
        GaussLegendre_weak_product_matrix!(P2,coord.ngrid,xreverse,wreverse,"P2",radau=true)

        Y00 = allocate_float(coord.ngrid, coord.ngrid, coord.ngrid)
        Y01 = allocate_float(coord.ngrid, coord.ngrid, coord.ngrid)
        Y10 = allocate_float(coord.ngrid, coord.ngrid, coord.ngrid)
        Y11 = allocate_float(coord.ngrid, coord.ngrid, coord.ngrid)
        Y20 = allocate_float(coord.ngrid, coord.ngrid, coord.ngrid)
        Y21 = allocate_float(coord.ngrid, coord.ngrid, coord.ngrid)
        Y30 = allocate_float(coord.ngrid, coord.ngrid, coord.ngrid)
        Y31 = allocate_float(coord.ngrid, coord.ngrid, coord.ngrid)
        GaussLegendre_weak_product_matrix!(Y00,coord.ngrid,xreverse,wreverse,"Y00",radau=true)
        GaussLegendre_weak_product_matrix!(Y01,coord.ngrid,xreverse,wreverse,"Y01",radau=true)
        GaussLegendre_weak_product_matrix!(Y10,coord.ngrid,xreverse,wreverse,"Y10",radau=true)
        GaussLegendre_weak_product_matrix!(Y11,coord.ngrid,xreverse,wreverse,"Y11",radau=true)
        GaussLegendre_weak_product_matrix!(Y20,coord.ngrid,xreverse,wreverse,"Y20",radau=true)
        GaussLegendre_weak_product_matrix!(Y21,coord.ngrid,xreverse,wreverse,"Y21",radau=true)
        GaussLegendre_weak_product_matrix!(Y30,coord.ngrid,xreverse,wreverse,"Y30",radau=true)
        GaussLegendre_weak_product_matrix!(Y31,coord.ngrid,xreverse,wreverse,"Y31",radau=true)
    else
        M2 = allocate_float(0, 0)
        S0 = allocate_float(0, 0)
        S1 = allocate_float(0, 0)
        K2 = allocate_float(0, 0)
        P1 = allocate_float(0, 0)
        P2 = allocate_float(0, 0)

        Y00 = allocate_float(0, 0, 0)
        Y01 = allocate_float(0, 0, 0)
        Y10 = allocate_float(0, 0, 0)
        Y11 = allocate_float(0, 0, 0)
        Y20 = allocate_float(0, 0, 0)
        Y21 = allocate_float(0, 0, 0)
        Y30 = allocate_float(0, 0, 0)
        Y31 = allocate_float(0, 0, 0)
    end
    return gausslegendre_base_info(false,Dmat, indefinite_integration_matrix, M0, M1, M2, S0,
                                   S1, K0, K1, K2, P0, P1, P2, D0, Y00, Y01, Y10, Y11,
                                   Y20, Y21, Y30, Y31)
end 

"""
A function that takes the first derivative in each element of `coord.grid`,
leaving the result (element-wise) in `coord.scratch_2d`.
"""
function elementwise_derivative!(coord, ff, gausslegendre::gausslegendre_info)
    df = coord.scratch_2d
    # define local variable nelement for convenience
    nelement = coord.nelement_local
    # check array bounds
    @debug_consistency_checks nelement == size(df,2) && coord.ngrid == size(df,1) || throw(BoundsError(df))
    
    # variable k will be used to avoid double counting of overlapping point
    k = 0
    j = 1 # the first element
    imin = coord.imin[j]-k
    # imax is the maximum index on the full grid for this (jth) element
    imax = coord.imax[j]        
    if coord.radau_first_element && coord.irank == 0 # differentiate this element with the Radau scheme
        @views mul!(df[:,j],gausslegendre.radau.Dmat[:,:],ff[imin:imax])
    else #differentiate using the Lobatto scheme
        @views mul!(df[:,j],gausslegendre.lobatto.Dmat[:,:],ff[imin:imax])
    end
    # transform back to the physical coordinate scale
    for i in 1:coord.ngrid
        df[i,j] /= coord.element_scale[j]
    end
    # calculate the derivative on each element
    @inbounds for j ∈ 2:nelement
        k = 1 
        imin = coord.imin[j]-k
        # imax is the maximum index on the full grid for this (jth) element
        imax = coord.imax[j]
        @views mul!(df[:,j],gausslegendre.lobatto.Dmat[:,:],ff[imin:imax])        
        # transform back to the physical coordinate scale
        for i in 1:coord.ngrid
            df[i,j] /= coord.element_scale[j]
        end
    end

    return nothing
end

"""
Wrapper function for element-wise derivatives with advection.
Note that Gauss-Legendre spectral the element method implemented here
does not use upwinding within an element.
"""
# Spectral element method does not use upwinding within an element
function elementwise_derivative!(coord, ff, speed, spectral::gausslegendre_info)
    return elementwise_derivative!(coord, ff, spectral)
end

"""
Function to perform interpolation on a single element.
"""
function single_element_interpolate!(result, newgrid, f, imin, imax, ielement, coord,
                                     gausslegendre::gausslegendre_base_info,
                                     derivative::Val{0})
    n_new = length(newgrid)

    i = 1
    other_nodes = @view coord.other_nodes[:,i,ielement]
    one_over_denominator = coord.one_over_denominator[i,ielement]
    this_f = f[i]
    for j ∈ 1:n_new
        result[j] = this_f * lagrange_poly_optimised(other_nodes, one_over_denominator, newgrid[j])
    end
    for i ∈ 2:coord.ngrid
        other_nodes = @view coord.other_nodes[:,i,ielement]
        one_over_denominator = coord.one_over_denominator[i,ielement]
        this_f = f[i]
        for j ∈ 1:n_new
            result[j] += this_f * lagrange_poly_optimised(other_nodes, one_over_denominator, newgrid[j])
        end
    end

    return nothing
end

"""
Function to carry out a 1D (global) mass matrix solve.
"""
# Evaluate first derivative of the interpolating function
function single_element_interpolate!(result, newgrid, f, imin, imax, ielement, coord,
                                     gausslegendre::gausslegendre_base_info,
                                     derivative::Val{1})
    n_new = length(newgrid)

    i = 1
    other_nodes = @view coord.other_nodes[:,i,ielement]
    one_over_denominator = coord.one_over_denominator[i,ielement]
    this_f = f[i]
    for j ∈ 1:n_new
        result[j] = this_f * lagrange_poly_derivative_optimised(other_nodes, one_over_denominator, newgrid[j])
    end
    for i ∈ 2:coord.ngrid
        other_nodes = @view coord.other_nodes[:,i,ielement]
        one_over_denominator = coord.one_over_denominator[i,ielement]
        this_f = f[i]
        for j ∈ 1:n_new
            result[j] += this_f * lagrange_poly_derivative_optimised(other_nodes, one_over_denominator, newgrid[j])
        end
    end

    return nothing
end

function fill_single_element_interpolation_matrix!(
             matrix_slice, newgrid, jelement, coord,
             gausslegendre::gausslegendre_base_info)
    n_new = length(newgrid)

    for j ∈ 1:coord.ngrid
        other_nodes = @view coord.other_nodes[:,j,jelement]
        one_over_denominator = coord.one_over_denominator[j,jelement]
        for i ∈ 1:n_new
            matrix_slice[i,j] = lagrange_poly_optimised(other_nodes, one_over_denominator, newgrid[i])
        end
    end

    return nothing
end

"""
Function to calculate the elemental anti-differentiation (integration) matrix.
This function forms the primitive
```math
F(x) = \\int^x_{x_{\\rm min}} f(x^\\prime) d x^\\prime
```
of the function 
```math
f(x) = \\sum_j f_j l_j(x),
```
where \$l_j(x)\$ is the \$j^{\\rm th}\$ Lagrange polynomial on the element and \$f_j = f(x_j)\$,
with \$x_j\$ \$j^{\\rm th}\$ collocation point on the element. We find \$F(x)\$ at the collocation
points on the element, giving a series of integrals to evaluate:
```math
F(x_i) = \\int^{x_i}_{-1} f(x^\\prime) d x^\\prime = \\sum_j f_j \\int^{x_i}_{-1} l_j(x^\\prime) d x^\\prime,
```
where we have used that \$x_{\\rm min} = -1\$ on the elemental grid.
Changing to a normalised coordinate \$y\$ suitable for Gaussian quadrature 
```math
x^\\prime = \\frac{x_i + 1}{2} y + \\frac{x_i - 1}{2}
```
we can write the operation in matrix form:
```math
F(x_i) = \\sum_{j}A_{ij}f_j,
```
with the matrix \$A_{ij}\$ defined by
```math
A_{ij} = \\left(\\frac{x_i + 1}{2}\\right) \\int^1_{-1} l_j \\left( \\frac{(x_i + 1)y + x_i - 1}{2} \\right) dy,
```
or in discretised form
```math
A_{ij} = \\left(\\frac{x_i + 1}{2}\\right) \\sum_k l_j \\left( \\frac{(x_i + 1)y_k + x_i - 1}{2} \\right) w_k,
```
with \$y_k\$ and \$w_k\$ Gauss-quadrature points and weights, respectively.
"""
function integration_matrix!(A::Array{Float64,2},x::Array{Float64,1},ngrid::Int64)
    nquad = 2*ngrid
    zz, wz = gausslegendre(nquad)
    # set lower limit
    for j in 1:ngrid
        A[1,j] = 0.0
    end
    # set the remaining points
    for i in 1:ngrid
        for j in 1:ngrid
            xi = x[i] # limit of ith integral
            scale = 0.5*(xi + 1)
            shift = 0.5*(xi - 1)
            # calculate the sum with Gaussian quadrature
            A[i,j] = 0.0
            for k in 1:nquad
                xval = scale*zz[k] + shift
                A[i,j] += scale*lagrange_poly(j,x,xval)*wz[k]
            end
        end
    end
    return nothing
end

"""
Formula for Gauss-Legendre-Lobatto differentiation matrix taken from p196 of Chpt `The Spectral Elemtent Method' of 
`Computational Seismology'. Heiner Igel First Edition. Published in 2017 by Oxford University Press.
Or https://doc.nektar.info/tutorials/latest/fundamentals/differentiation/fundamentals-differentiationch2.html

    D -- differentiation matrix 
    x -- Gauss-Legendre-Lobatto points in [-1,1]
    ngrid -- number of points per element (incl. boundary points)

Note that D has does not include a scaling factor
"""
function gausslobattolegendre_differentiation_matrix!(D::Array{mk_float,2},x::Array{mk_float,1},ngrid::mk_int)
    D[:,:] .= 0.0
    for ix in 1:ngrid
        for ixp in 1:ngrid
            if !(ix == ixp)
                D[ix,ixp] = (Pl(x[ix],ngrid-1)/Pl(x[ixp],ngrid-1))/(x[ix]-x[ixp])
            end
        end
    end
    # uncomment for analytical diagonal values 
    #D[1,1] = -0.25*(ngrid - 1)*ngrid
    #D[ngrid,ngrid] = 0.25*(ngrid - 1)*ngrid
    #for ix in 1:ngrid-1
    #   D[ix,ix] = 0.0
    #end
    # get diagonal values from sum of nonzero off diagonal values 
    for ix in 1:ngrid
        D[ix,ix] = -sum(D[ix,:])
    end 
    return nothing
end

"""
Formula for Gauss-Legendre-Radau differentiation matrix taken from
https://doc.nektar.info/tutorials/latest/fundamentals/differentiation/fundamentals-differentiationch2.html

    D -- differentiation matrix 
    x -- Gauss-Legendre-Radau points in [-1,1)
    ngrid -- number of points per element (incl. boundary points)

Note that D has does not include a scaling factor
"""
function gaussradaulegendre_differentiation_matrix!(D::Array{mk_float,2},x::Array{mk_float,1},ngrid::Int64)
    D[:,:] .= 0.0
    for ix in 1:ngrid
        for ixp in 1:ngrid
            if !(ix == ixp)
                D[ix,ixp] = (Pl(x[ix],ngrid-1)/Pl(x[ixp],ngrid-1))*((1.0 - x[ixp])/(1.0 - x[ix]))/(x[ix]-x[ixp])
            end
        end
    end
    # uncomment for analytical diagonal values 
    #D[1,1] = -0.25*(ngrid - 1)*(ngrid + 1)
    #for ix in 2:ngrid
    #   D[ix,ix] = 0.5/(1.0 - x[ix])
    #end
    # get diagonal values from sum of nonzero off diagonal values 
    for ix in 1:ngrid
        D[ix,ix] = -sum(D[ix,:])
    end
    
    # get into correct order for a grid on (-1,1]
    Dreverse = copy(D)
    for ix in 1:ngrid
        for ixp in 1:ngrid
            Dreverse[ngrid-ix+1,ngrid-ixp+1] = -D[ix,ixp]
        end
    end
    D .= Dreverse
    return nothing
end

"""
Gauss-Legendre derivative at arbitrary x values, for boundary condition on Radau points.
    D0 -- the vector
    xj -- the x location where the derivative is evaluated 
    ngrid -- number of points in x
    x -- the grid from -1, 1
Note that D0 is not scaled to the physical grid with a scaling factor.
"""
function GaussLegendre_derivative_vector!(D0,xj,ngrid,x,wgts;radau=false)
    # coefficient in expansion of 
    # lagrange polys in terms of Legendre polys
    gamma = allocate_float(ngrid)
    for i in 1:ngrid-1
        gamma[i] = Legendre_h_n(i-1)
    end
    if radau
        gamma[ngrid] = Legendre_h_n(ngrid-1)
    else
        gamma[ngrid] = 2.0/(ngrid - 1)
    end
    
    @. D0 = 0.0
    for i in 1:ngrid
        for k in 1:ngrid
            D0[i] += wgts[i]*Pl(x[i],k-1)*dnPl(xj,k-1,1)/gamma[k]
        end
    end
    # set `diagonal' value
    D0[1] = 0.0
    D0[1] = -sum(D0[:])
    #@. D0 *= 2.0*float(nelement_global)/L
end

"""
Result of the inner product of Legendre polynomials of order k.
"""
function Legendre_h_n(k)
    h_n = 2.0/(2.0*k + 1)
    return h_n
end 


"""
Assign abitrary weak inner product matrix `Q` on a 1D line with Jacobian equal to 1
matrix `Q` acts on a single vector `x` such that `y = Q * x` is also a vector.

We use a projection onto Gauss-Legendre polynomials to carry out the calculation
in two steps (see, e.g, S. A. Teukolsky, Short note on the mass matrix for Gauss–Lobatto grid
points, J. Comput. Phys. 283 (2015) 408–413. https://doi.org/10.1016/j.jcp.2014.12.012).
First, we write the desired matrix elements in terms of Legendre polynomials
```math
   l_i(x) = \\sum_j \\frac{P_j(x)P_j(x_i)w_i}{\\gamma_j}
```
with \$w_i\$ the weights from an integration on the Gauss-Legendre-Lobatto (or Radau) points \$x_i\$,
i.e., 
```math 
 \\int^1_{-1} f(x) d x = \\sum_{i} f(x_i)w_i,
```
and \$\\gamma_j = \\sum_k w_k P_j(x_k)P_j(x_k)\$ the numerical inner-product.
Then, a matrix element can be expressed in integrals over Legendre polynomials
rather than Lagrange polynomials, i.e.,
```math
   M_{ij} = \\int^1_{-1} l_i(x)l_j(x) d x = \\sum_{mn} \\frac{w_m P_m(x_i) w_n P_n(x_j)}{\\gamma_m\\gamma_n} \\int^{1}_{-1} P_m(x)P_n(x) d x. 
```
Defining 
```math
  A_{mn} = \\int^{1}_{-1} P_m(x)P_n(x) d x, 
```
we can thus write
```math
   M_{ij} = \\sum_{mn} \\frac{w_m P_m(x_i) w_n P_n(x_j)}{\\gamma_m\\gamma_n} A_{mn}. 
```
We can use a quadrature which yields exact results (to machine precision)
to evaluate \$A_{mn}\$ using fast library functions for the Legendre polynomials,
and then carry out the sum \$\\sum_{mn}\$ to obtain exact results (to machine-precision).
Here we use a Gauss-Legendre integration quadrature with exact results up to 
polynomials with order \$k_{max} = 4N +1\$, with \$N=\$`ngrid` and the highest order polynomial product
that we integrate is \$P_{N-1}(x)P_{N-1}(x)x^2\$, which has order \$k=2N < k_{max}\$.

"""
function GaussLegendre_weak_product_matrix!(QQ::Array{mk_float,2},ngrid,x,wgts,option;radau=false)
    # coefficient in expansion of 
    # lagrange polys in terms of Legendre polys
    gamma = allocate_float(ngrid)
    for i in 1:ngrid-1
        gamma[i] = Legendre_h_n(i-1)
    end
    if radau
        gamma[ngrid] = Legendre_h_n(ngrid-1)
    else
        gamma[ngrid] = 2.0/(ngrid - 1)
    end
    # appropriate inner product of Legendre polys
    # definition depends on required matrix 
    # for M0: AA = < P_i P_j >
    # for M1: AA = < P_i P_j x >
    # for M2: AA = < P_i P_j x^2 >
    # for S0: AA = -< P'_i P_j >
    # for S1: AA = -< P'_i P_j x >
    # for K0: AA = -< P'_i P'_j >
    # for K1: AA = -< P'_i P'_j x >
    # for K2: AA = -< P'_i P'_j x^2 >
    # for P0: AA = < P_i P'_j >
    # for P1: AA = < P_i P'_j x >
    # for P2: AA = < P_i P'_j x^2 >
    AA = allocate_float(ngrid, ngrid)
    nquad = 2*ngrid
    zz, wz = gausslegendre(nquad)
    @. AA = 0.0
    if option == "M0"
        for j in 1:ngrid
            for i in 1:ngrid
                for k in 1:nquad
                    AA[i,j] += wz[k]*Pl(zz[k],i-1)*Pl(zz[k],j-1)
                end
            end
        end
    elseif option == "M1"
        for j in 1:ngrid
            for i in 1:ngrid
                for k in 1:nquad
                    AA[i,j] += zz[k]*wz[k]*Pl(zz[k],i-1)*Pl(zz[k],j-1)
                end
            end
        end
    elseif option == "M2"
        for j in 1:ngrid
            for i in 1:ngrid
                for k in 1:nquad
                    AA[i,j] += (zz[k]^2)*wz[k]*Pl(zz[k],i-1)*Pl(zz[k],j-1)
                end
            end
        end
    elseif option == "S0"
        for j in 1:ngrid
            for i in 1:ngrid
                for k in 1:nquad
                    AA[i,j] -= wz[k]*dnPl(zz[k],i-1,1)*Pl(zz[k],j-1)
                end
            end
        end
    elseif option == "S1"
        for j in 1:ngrid
            for i in 1:ngrid
                for k in 1:nquad
                    AA[i,j] -= zz[k]*wz[k]*dnPl(zz[k],i-1,1)*Pl(zz[k],j-1)
                end
            end
        end
    elseif option == "K0"
        for j in 1:ngrid
            for i in 1:ngrid
                for k in 1:nquad
                    AA[i,j] -= wz[k]*dnPl(zz[k],i-1,1)*dnPl(zz[k],j-1,1)
                end
            end
        end
    elseif option == "K1"
        for j in 1:ngrid
            for i in 1:ngrid
                for k in 1:nquad
                    AA[i,j] -= zz[k]*wz[k]*dnPl(zz[k],i-1,1)*dnPl(zz[k],j-1,1)
                end
            end
        end
    elseif option == "K2"
        for j in 1:ngrid
            for i in 1:ngrid
                for k in 1:nquad
                    AA[i,j] -= (zz[k]^2)*wz[k]*dnPl(zz[k],i-1,1)*dnPl(zz[k],j-1,1)
                end
            end
        end
    elseif option == "P0"
        for j in 1:ngrid
            for i in 1:ngrid
                for k in 1:nquad
                    AA[i,j] += wz[k]*Pl(zz[k],i-1)*dnPl(zz[k],j-1,1)
                end
            end
        end
    elseif option == "P1"
        for j in 1:ngrid
            for i in 1:ngrid
                for k in 1:nquad
                    AA[i,j] += zz[k]*wz[k]*Pl(zz[k],i-1)*dnPl(zz[k],j-1,1)
                end
            end
        end
    elseif option == "P2"
        for j in 1:ngrid
            for i in 1:ngrid
                for k in 1:nquad
                    AA[i,j] += (zz[k]^2)*wz[k]*Pl(zz[k],i-1)*dnPl(zz[k],j-1,1)
                end
            end
        end
    end
    
    QQ .= 0.0
    for j in 1:ngrid
        for i in 1:ngrid
            for l in 1:ngrid
                for k in 1:ngrid
                    QQ[i,j] += wgts[i]*wgts[j]*Pl(x[i],k-1)*Pl(x[j],l-1)*AA[k,l]/(gamma[k]*gamma[l])
                end
            end
        end
    end
    return nothing
end

"""
Assign abitrary weak (nonlinear) inner product matrix `Q` on a 1D line with Jacobian equal to 1.
matrix `Q` acts on two vectors `x1` and `x2` such that the quadratic form 
`y = x1 * Q * x2` is also a vector. See documentation of corresponding function
for linear inner product matrices.
"""
function GaussLegendre_weak_product_matrix!(QQ::Array{mk_float,3},ngrid,x,wgts,option;radau=false)
    # coefficient in expansion of 
    # lagrange polys in terms of Legendre polys
    gamma = allocate_float(ngrid)
    for i in 1:ngrid-1
        gamma[i] = Legendre_h_n(i-1)
    end
    if radau
        gamma[ngrid] = Legendre_h_n(ngrid-1)
    else
        gamma[ngrid] = 2.0/(ngrid - 1)
    end
    # appropriate inner product of Legendre polys
    # definition depends on required matrix 
    # for Y00: AA = < P_i P_j P_k >
    # for Y01: AA = < P_i P_j P_k x >
    # for Y10: AA = < P_i P_j P'_k >
    # for Y11: AA = < P_i P_j P'_k x >
    # for Y20: AA = < P_i P'_j P'_k >
    # for Y21: AA = < P_i P'_j P'_k x >
    # for Y31: AA = < P_i P'_j P_k x >
    # for Y30: AA = < P_i P'_j P_k >
    AA = allocate_float(ngrid, ngrid, ngrid)
    nquad = 2*ngrid
    zz, wz = gausslegendre(nquad)
    @. AA = 0.0
    if option == "Y00"
        for k in 1:ngrid
            for j in 1:ngrid
                for i in 1:ngrid
                    for q in 1:nquad
                        AA[i,j,k] += wz[q]*Pl(zz[q],i-1)*Pl(zz[q],j-1)*Pl(zz[q],k-1)
                    end
                end
            end
        end
    elseif option == "Y01"
        for k in 1:ngrid
            for j in 1:ngrid
                for i in 1:ngrid
                    for q in 1:nquad
                        AA[i,j,k] += zz[q]*wz[q]*Pl(zz[q],i-1)*Pl(zz[q],j-1)*Pl(zz[q],k-1)
                    end
                end
            end
        end
    elseif option == "Y10"
        for k in 1:ngrid
            for j in 1:ngrid
                for i in 1:ngrid
                    for q in 1:nquad
                        AA[i,j,k] += wz[q]*Pl(zz[q],i-1)*Pl(zz[q],j-1)*dnPl(zz[q],k-1,1)
                    end
                end
            end
        end
    elseif option == "Y11"
        for k in 1:ngrid
            for j in 1:ngrid
                for i in 1:ngrid
                    for q in 1:nquad
                        AA[i,j,k] += zz[q]*wz[q]*Pl(zz[q],i-1)*Pl(zz[q],j-1)*dnPl(zz[q],k-1,1)
                    end
                end
            end
        end
    elseif option == "Y20"
        for k in 1:ngrid
            for j in 1:ngrid
                for i in 1:ngrid
                    for q in 1:nquad
                        AA[i,j,k] += wz[q]*Pl(zz[q],i-1)*dnPl(zz[q],j-1,1)*dnPl(zz[q],k-1,1)
                    end
                end
            end
        end
    elseif option == "Y21"
        for k in 1:ngrid
            for j in 1:ngrid
                for i in 1:ngrid
                    for q in 1:nquad
                        AA[i,j,k] += zz[q]*wz[q]*Pl(zz[q],i-1)*dnPl(zz[q],j-1,1)*dnPl(zz[q],k-1,1)
                    end
                end
            end
        end
    elseif option == "Y31"
        for k in 1:ngrid
            for j in 1:ngrid
                for i in 1:ngrid
                    for q in 1:nquad
                        AA[i,j,k] += zz[q]*wz[q]*Pl(zz[q],i-1)*dnPl(zz[q],j-1,1)*Pl(zz[q],k-1)
                    end
                end
            end
        end
    elseif option == "Y30"
        for k in 1:ngrid
            for j in 1:ngrid
                for i in 1:ngrid
                    for q in 1:nquad
                        AA[i,j,k] += wz[q]*Pl(zz[q],i-1)*dnPl(zz[q],j-1,1)*Pl(zz[q],k-1)
                    end
                end
            end
        end
    end
    
    QQ .= 0.0
    for k in 1:ngrid
        for j in 1:ngrid
            for i in 1:ngrid
                for l in 1:ngrid
                    for m in 1:ngrid
                        for n in 1:ngrid
                            QQ[i,j,k] += wgts[i]*wgts[j]*wgts[k]*Pl(x[i],n-1)*Pl(x[j],m-1)*Pl(x[k],l-1)*AA[n,m,l]/(gamma[n]*gamma[m]*gamma[l])
                        end
                    end
                end
            end
        end
    end
    return nothing
end

"""
Function for computing the scale factor on a grid with uniformed spaced element boundaries.
Unused.
"""
function scale_factor_func(L,nelement_global)
    return 0.5*L/float(nelement_global)
end

"""
Function for computing the shift factor on a grid with uniformed spaced element boundaries.
Unused.
"""
function shift_factor_func(L,nelement_global,nelement_local,irank,ielement_local)
    #ielement_global = ielement_local # for testing + irank*nelement_local
    ielement_global = ielement_local + irank*nelement_local # proper line for future distributed memory MPI use
    shift = L*((float(ielement_global)-0.5)/float(nelement_global) - 0.5)
    return shift
end

"""
Function for finding the elemental index in the global distributed-memory grid.
Distributed-memory for global finite-element operators is not yet supported.
"""
function ielement_global_func(nelement_local,irank,ielement_local)
    return ielement_global = ielement_local + irank*nelement_local
end

"""
Function for setting up the full Gauss-Legendre-Lobatto
grid and collocation point weights.
"""
function scaled_gauss_legendre_lobatto_grid(ngrid, nelement_local, n_local, element_scale, element_shift, imin, imax)
    # get Gauss-Legendre-Lobatto points and weights on [-1,1]
    x, w = gausslobatto(ngrid)
    # grid and weights arrays
    grid = allocate_float(n_local)
    wgts = allocate_float(n_local)
    wgts .= 0.0
    #integer to deal with the overlap of element boundaries
    k = 1
    @inbounds for j in 1:nelement_local
        # element_scale[j]
        # element_shift[j]
        # factor with maps [-1,1] -> a subset of [-L/2, L/2]
        @. grid[imin[j]:imax[j]] = element_scale[j]*x[k:ngrid] + element_shift[j]
        
        # calculate the weights
        # remembering on boundary points to include weights
        # from both left and right elements
        #println(imin[j]," ",imax[j])
        @. wgts[imin[j] - k + 1:imax[j]] += element_scale[j]*w[1:ngrid] 
        
        k = 2        
    end
    return grid, wgts
end

"""
Function for setting up the full Gauss-Legendre-Radau
grid and collocation point weights.
"""
function scaled_gauss_legendre_radau_grid(ngrid, nelement_local, n_local, element_scale, element_shift, imin, imax, irank)
    # get Gauss-Legendre-Lobatto points and weights on [-1,1]
    x_lob, w_lob = gausslobatto(ngrid)
    # get Gauss-Legendre-Radau points and weights on [-1,1)
    x_rad, w_rad = gaussradau(ngrid)
    # transform to a Gauss-Legendre-Radau grid on (-1,1]
    x_rad, w_rad = -reverse(x_rad), reverse(w_rad)#
    # grid and weights arrays
    grid = allocate_float(n_local)
    wgts = allocate_float(n_local)
    wgts .= 0.0
    if irank == 0
        # for 1st element, fill in with Gauss-Legendre-Radau points
        j = 1
        # element_scale[j]
        # element_shift[j]
        # factor with maps [-1,1] -> a subset of [-L/2, L/2]
        @. grid[imin[j]:imax[j]] = element_scale[j]*x_rad[1:ngrid] + element_shift[j]
        @. wgts[imin[j]:imax[j]] += element_scale[j]*w_rad[1:ngrid]       
        #integer to deal with the overlap of element boundaries
        k = 2
        @inbounds for j in 2:nelement_local
            # element_scale[j]
            # element_shift[j]
            # factor with maps [-1,1] -> a subset of [-L/2, L/2]
            @. grid[imin[j]:imax[j]] = element_scale[j]*x_lob[k:ngrid] + element_shift[j]
            @. wgts[imin[j] - k + 1:imax[j]] += element_scale[j]*w_lob[1:ngrid]         
        end
    else # all elements are Gauss-Legendre-Lobatto
        #integer to deal with the overlap of element boundaries
        k = 1
        @inbounds for j in 1:nelement_local
            # element_scale[j]
            # element_shift[j]
            # factor with maps [-1,1] -> a subset of [-L/2, L/2]
            @. grid[imin[j]:imax[j]] = element_scale[j]*x_lob[k:ngrid] + element_shift[j]
            @. wgts[imin[j] - k + 1:imax[j]] += element_scale[j]*w_lob[1:ngrid]            
            k = 2 
        end
    end
    return grid, wgts
end

end
