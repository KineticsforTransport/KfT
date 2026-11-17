"""
Lagrange polynomials can be useful for finite element methods on any set of basis points,
as they give a representation of the interpolating function within an element whose
coefficients are just the function values at the grid points.

This module collects some functions related to the use of Lagrange polynomials, to avoid
their being scattered (and possibly duplicated) in other modules.
"""
module lagrange_polynomials

export lagrange_poly, lagrange_poly_optimised, lagrange_poly_derivative_optimised,
       lagrange_weak_product_matrix

import FastGaussQuadrature

using ..array_allocation: allocate_float
using ..type_definitions: mk_float

"""
Lagrange polynomial
args:
j - index of l_j from list of nodes
x_nodes - array of x node values
x - point where interpolated value is returned
"""
function lagrange_poly(j,x_nodes,x)
    # get number of nodes
    n = size(x_nodes,1)
    # location where l(x0) = 1
    x0 = x_nodes[j]
    # evaluate polynomial
    poly = 1.0
    for i in 1:j-1
            poly *= (x - x_nodes[i])/(x0 - x_nodes[i])
    end
    for i in j+1:n
            poly *= (x - x_nodes[i])/(x0 - x_nodes[i])
    end
    return poly
end

"""
    lagrange_poly_optimised(other_nodes, one_over_denominator, x)

Optimised version of Lagrange polynomial calculation, making use of pre-calculated quantities.

`other_nodes` is a vector of the grid points in this element where this Lagrange
polynomial is zero (the other nodes than the one where it is 1).

`one_over_denominator` is `1/prod(x0 - n for n ∈ other_nodes)` where `x0` is the grid
point where this Lagrange polynomial is 1.

`x` is the point to evaluate the Lagrange polynomial at.
"""
function lagrange_poly_optimised(other_nodes, one_over_denominator, x)
    return prod(x - n for n ∈ other_nodes) * one_over_denominator
end

"""
    lagrange_poly_derivative_optimised(other_nodes, one_over_denominator, x)

Optimised calculation of the first derivative of a Lagrange polynomial, making use of
pre-calculated quantities.

`other_nodes` is a vector of the grid points in this element where this Lagrange
polynomial is zero (the other nodes than the one where it is 1).

`one_over_denominator` is `1/prod(x0 - n for n ∈ other_nodes)` where `x0` is the grid
point where this Lagrange polynomial is 1.

`x` is the point to evaluate the Lagrange polynomial at.
"""
function lagrange_poly_derivative_optimised(other_nodes, one_over_denominator, x)
    result = 0.0
    k = length(other_nodes)
    # Is there a more efficient way of doing this? Not a big deal for now because this
    # function will only be used to calculate a preconditioner matrix, which is done
    # rarely.
    for i ∈ 1:k
        result += prod(x - other_nodes[j] for j ∈ 1:k if j ≠ i)
    end
    result *= one_over_denominator
    return result
end

"""
    lagrange_weak_product_matrix!(QQ::Array{mk_float,2}, ngrid,
                                  lagrange_nodes::AbstractVector{mk_float},
                                  lagrange_other_nodes::AbstractMatrix{mk_float},
                                  lagrande_one_over_denominator::AbstractVector{mk_float},
                                  option; radau=false)

Assign abitrary weak inner product matrix `Q` on a 1D line with Jacobian equal to 1
matrix `Q` acts on a single vector `x` such that `y = Q * x` is also a vector.

We evaluate the degree-(ngrid-1) Lagrange polynomials on a Gaussian quadrature grid with
ngrid points that can exactly integrate polynomials up to degree (2*ngrid-1), which is
enough here as the maximum polynomial degree is 2*(ngrid-1)+1 (in M1).

Note that the `lagrange_nodes` here should be defined on the range [-1,1], as this 'local
matrix' will be scaled later to fit different elements. The `lagrange_other_nodes` and
`lagrange_one_over_denominator` should have corresponding values.
"""
function lagrange_weak_product_matrix!(QQ::Array{mk_float,2}, ngrid,
                                       lagrange_nodes::AbstractVector{mk_float},
                                       lagrange_other_nodes::AbstractMatrix{mk_float},
                                       lagrange_one_over_denominator::AbstractVector{mk_float},
                                       option)
    nquad = ngrid
    quadgrid, quadweights = FastGaussQuadrature.gausslegendre(nquad)

    if option == "M0"
        # 'Mass matrix', ∫ L_i(x) L_j(x) dx

        # Evaluate Lagrange polynomials corresponding to `lagrange_nodes` on quadrature
        # grid.
        poly_values = allocate_float(nquad, ngrid)
        for igrid ∈ 1:ngrid, jquad ∈ 1:nquad
            poly_values[jquad,igrid] =
                @views lagrange_poly_optimised(lagrange_other_nodes[:,igrid],
                                               lagrange_one_over_denominator[igrid],
                                               quadgrid[jquad])
        end

        # Evaluate the integrals of pairs of polynomials using Gaussian quadrature.
        for j ∈ 1:ngrid, i ∈ 1:j
            # q = ∑_k L_i(x_k) L_j(x_k) w_k
            # where x_k, w_k are the values in quadgrid and quadweights.
            QQ[i,j] = @views mapreduce(*, +, poly_values[:,i], poly_values[:,j], quadweights)
            if i != j
                # The mass matrix is symmetric, so just copy (i,j) into (j,i)
                QQ[j,i] = QQ[i,j]
            end
        end
    elseif option == "M1"
        # 'Mass matrix type 1', ∫ x L_i(x) L_j(x) dx

        # Evaluate Lagrange polynomials corresponding to `lagrange_nodes` on quadrature
        # grid.
        poly_values = allocate_float(nquad, ngrid)
        for igrid ∈ 1:ngrid, jquad ∈ 1:nquad
            poly_values[jquad,igrid] =
                @views lagrange_poly_optimised(lagrange_other_nodes[:,igrid],
                                               lagrange_one_over_denominator[igrid],
                                               quadgrid[jquad])
        end

        # Evaluate the integrals of pairs of polynomials using Gaussian quadrature.
        for j ∈ 1:ngrid, i ∈ 1:j
            # q = ∑_k L_i(x_k) L_j(x_k) w_k
            # where x_k, w_k are the values in quadgrid and quadweights.
            QQ[i,j] = @views mapreduce(*, +, quadgrid, poly_values[:,i], poly_values[:,j],
                                       quadweights)
            if i != j
                # The mass matrix is symmetric, so just copy (i,j) into (j,i)
                QQ[j,i] = QQ[i,j]
            end
        end
    elseif option == "K0"
        # 'Stiffness matrix', ∫ ∂L_i(x)/∂x ∂L_j(x)/∂x dx

        # Evaluate Lagrange polynomials derivatives corresponding to `lagrange_nodes` on
        # quadrature grid.
        poly_derivative_values = allocate_float(nquad, ngrid)
        for igrid ∈ 1:ngrid, jquad ∈ 1:nquad
            poly_derivative_values[jquad,igrid] =
                @views lagrange_poly_derivative_optimised(lagrange_other_nodes[:,igrid],
                                                          lagrange_one_over_denominator[igrid],
                                                          quadgrid[jquad])
        end

        # Evaluate the integrals of pairs of polynomials using Gaussian quadrature.
        for j ∈ 1:ngrid, i ∈ 1:j
            # q = ∑_k L_i(x_k) L_j(x_k) w_k
            # where x_k, w_k are the values in quadgrid and quadweights.
            QQ[i,j] = - @views mapreduce(*, +, poly_derivative_values[:,i],
                                         poly_derivative_values[:,j], quadweights)
            if i != j
                # The mass matrix is symmetric, so just copy (i,j) into (j,i)
                QQ[j,i] = QQ[i,j]
            end
        end
    elseif option == "K1"
        # 'Stiffness matrix type 1', ∫ x ∂L_i(x)/∂x ∂L_j(x)/∂x dx

        # Evaluate Lagrange polynomials derivatives corresponding to `lagrange_nodes` on
        # quadrature grid.
        poly_derivative_values = allocate_float(nquad, ngrid)
        for igrid ∈ 1:ngrid, jquad ∈ 1:nquad
            poly_derivative_values[jquad,igrid] =
                @views lagrange_poly_derivative_optimised(lagrange_other_nodes[:,igrid],
                                                          lagrange_one_over_denominator[igrid],
                                                          quadgrid[jquad])
        end

        # Evaluate the integrals of pairs of polynomials using Gaussian quadrature.
        for j ∈ 1:ngrid, i ∈ 1:j
            # q = ∑_k L_i(x_k) L_j(x_k) w_k
            # where x_k, w_k are the values in quadgrid and quadweights.
            QQ[i,j] = - @views mapreduce(*, +, quadgrid, poly_derivative_values[:,i],
                                         poly_derivative_values[:,j], quadweights)
            if i != j
                # The mass matrix is symmetric, so just copy (i,j) into (j,i)
                QQ[j,i] = QQ[i,j]
            end
        end
    elseif option == "P0"
        # 'Mixed matrix', ∫ L_i(x) ∂L_j(x)/∂x dx

        # Evaluate Lagrange polynomials corresponding to `lagrange_nodes` on quadrature
        # grid.
        poly_values = allocate_float(nquad, ngrid)
        for igrid ∈ 1:ngrid, jquad ∈ 1:nquad
            poly_values[jquad,igrid] =
                @views lagrange_poly_optimised(lagrange_other_nodes[:,igrid],
                                               lagrange_one_over_denominator[igrid],
                                               quadgrid[jquad])
        end

        # Evaluate Lagrange polynomials derivatives corresponding to `lagrange_nodes` on
        # quadrature grid.
        poly_derivative_values = allocate_float(nquad, ngrid)
        for igrid ∈ 1:ngrid, jquad ∈ 1:nquad
            poly_derivative_values[jquad,igrid] =
                @views lagrange_poly_derivative_optimised(lagrange_other_nodes[:,igrid],
                                                          lagrange_one_over_denominator[igrid],
                                                          quadgrid[jquad])
        end

        # Evaluate the integrals of pairs of polynomials using Gaussian quadrature.
        for j ∈ 1:ngrid, i ∈ 1:ngrid
            # q = ∑_k L_i(x_k) L_j(x_k) w_k
            # where x_k, w_k are the values in quadgrid and quadweights.
            QQ[i,j] = - @views mapreduce(*, +, poly_values[:,i],
                                         poly_derivative_values[:,j], quadweights)
        end
    else
        error("Unsupported option=\"$option\".")
    end
end

end
