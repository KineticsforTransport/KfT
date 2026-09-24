using CommonSolve: solve
using CurveFit
using makie_post_processing
using makie_post_processing.CairoMakie: Figure, Axis, Legend, save, lines!, scatter!
using makie_post_processing.LaTeXStrings
using SpecialFunctions: dawson

# The following copied/modified from moment_kinetics/test/harrisonthompson.jl:
# Analytic solution given by implicit equation
#   z = ± C * D(sqrt(-phi))
# with +'ve for z>0 and -'ve for z<0, where D() is the 'Dawson function'
#   D(x) = exp(-x^2) ∫_0^x exp(y^2) dy
#
# Note the derivative of the dawson function is
Dprime(x) = -2.0 * x * dawson(x) + 1.0

# Need to find the 'magic' constant C for which dphi/dz=±∞ at z=±0.5.
# Equivalently,
#   0 = dz/dphi
#     = ± C * (-1) / (2 * sqrt(-phi)) * D'(sqrt(phi))
# where
#  1/2 = C * D(sqrt(-phi))
# the constants are non-zero, so must be D'(sqrt(phi))=0,
#   0 = -2 * sqrt(-phi) * D(sqrt(-phi)) + 1
#     = -2 * sqrt(-phi) * (1/2C) + 1
#   C = sqrt(-phi(±1/2))
# So evaluating at z=1/2, defining phib=phi(1/2)
#   1/2 = sqrt(-phib) * D(sqrt(-phib))
#   1/2 = C * D(C)
# We can determine the corresponding value of C by fixed-point iteration.
function get_magic_C()
    C = 1.0
    old_C = -1.0
    update_func(x) = 0.5 / dawson(x)
    while abs(C - old_C) > 1.0e-14
        old_C = C
        C = update_func(C)
    end
    return C
end

function findphi_bisection(z, C)
    z = abs(z)
    max_X = C
    min_X = 0.0
    X = nothing
    while (max_X - min_X) > 1.0e-13
        mid_X = 0.5 * (min_X + max_X)
        testval = C * dawson(mid_X)
        if testval == z
            X = mid_X
            break
        elseif testval < z
            min_X = mid_X
        else
            max_X = mid_X
        end
    end
    if X === nothing
        X = 0.5 * (min_X + max_X)
    end
    return - X^2
end

# Note, don't use Newton iteration as it is unstable near the sheath entrances.
maxits = 1000000

function newton(f, fprime, z, args...)
    phi = -1.0e-6
    count = 0
    while abs(f(phi, z, args...)) > 1.e-14
        delta_phi = - f(phi, z, args...) / fprime(phi, z, args...)
        if phi + delta_phi > 0.0
            phi = phi + 0.01 * delta_phi
        else
            phi = phi + delta_phi
        end
        if phi > -eps()
            # phi must be negative, but iteration might overshoot
            phi = -eps()
        end
        count += 1
        if count > maxits
            error("Reached maximum iterations for z=$z with phi=$phi")
        end
    end
    return phi
end

# want to find phi such that f(phi) is zero
fnegative(phi, z, C) = - z - C * dawson(sqrt(-phi))
fpositive(phi, z, C) = - z + C * dawson(sqrt(-phi))
# derivative of f for Newton iteration
fprimenegative(phi, z, C) = C / 2 / sqrt(-phi) * Dprime(sqrt(-phi))
fprimepositive(phi, z, C) = - C / 2 / sqrt(-phi) * Dprime(sqrt(-phi))
function findphi(z, C)
    if abs(z) == 0.5
        # Boundary point - iteration would fail but we know the value would be -C^2
        return -C^2
    elseif z < - eps()
        return newton(fnegative, fprimenegative, z, C)
    elseif z > eps()
        return newton(fpositive, fprimepositive, z, C)
    else
        return 0.0
    end
end

function HarrisonThompson_comparison(run_names...)
    run_info = get_run_info(run_names...)
    dir_name = "comparison_plots"

    C = get_magic_C()

    fig = Figure(; size=(1200,800))
    ax = Axis(fig[1,1]; xlabel=L"$z$", ylabel=L"$\phi$")

    near_wall_fig = Figure(; size=(1200,800))
    near_wall_ax = Axis(near_wall_fig[1,1]; xlabel=L"$\sqrt{z - z_\mathrm{wall}}$", ylabel=L"$\phi$")

    HT_gridnum = 10000
    HT_z = (-HT_gridnum:HT_gridnum) ./ (2 * HT_gridnum)
    HT_phi = findphi_bisection.(HT_z, C)

    # Offset so phi is zero at the (lower) sheath entrance for all plots.
    HT_phi .-= HT_phi[1]

    lines!(ax, HT_z, HT_phi; color=:black, label="Harrison-Thompson analytic")

    filter_inds = HT_z .< -0.49
    sqrt_delta_HT_z = sqrt.(HT_z[filter_inds] .+ 0.5)
    lines!(near_wall_ax, sqrt_delta_HT_z, HT_phi[filter_inds]; color=:black, label="Harrison-Thompson analytic")

    # Plot phi profiles from simulations.
    for ri ∈ run_info
        data = get_variable(ri, "phi"; it=-1, ir=1)
        z = ri.z.grid

        # Offset so phi is zero at the (lower) sheath entrance for all plots.
        data .-= data[1]

        plot_vs_z(ri, "phi"; ax, data, linestyle=:dash, label=ri.run_name)

        filter_inds = z .< -0.49
        sqrt_delta_z = sqrt.(z[filter_inds] .+ 0.5)
        filtered_data = data[filter_inds]
        s = scatter!(near_wall_ax, sqrt_delta_z, filtered_data; label=ri.run_name)

        # Fit straight line to second half of data, to visualise deviation from sqrt(z)
        # behaviour near sheath entrance.
        n = length(sqrt_delta_z)
        x = sqrt_delta_z[n÷2:end]
        y = filtered_data[n÷2:end]
        prob = CurveFitProblem(x, y)
        sol = solve(prob, LinearCurveFitAlgorithm())
        plotx = collect(0.0:0.01:0.1)
        #lines!(near_wall_ax, plotx, sol(plotx); linestyle=:dot, color=s.color)
lines!(near_wall_ax, plotx, sol(plotx); linestyle=:dot, color=s.color)
println("plotx=$plotx, sol=", sol(plotx), " x=$x, y=$y n=$n")
    end

    Legend(fig[2,1], ax; tellwidth=false, tellheight=true)
    save(joinpath(dir_name, "HarrisonThompson.pdf"), fig)

    Legend(near_wall_fig[2,1], near_wall_ax; tellwidth=false, tellheight=true)
    save(joinpath(dir_name, "HarrisonThompson-near-wall.pdf"), near_wall_fig)

    return fig, ax, near_wall_fig, near_wall_ax
end
