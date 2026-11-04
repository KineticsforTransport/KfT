"""
Nonlinear solvers, using Jacobian-free Newton-Krylov methods.

These solvers use an outer Newton iteration. Each step of the Newton iteration requires a
linear solve of the Jacobian. An 'inexact Jacobian' method is used, and the GMRES method
(GMRES is a type of Krylov solver) is used to (approximately) solve the (approximate)
linear system.

!!! warning "parallelisation"
    This module uses shared- and distributed-memory parallelism, so the functions in it
    should not be called inside any kind of parallelised loop. This restriction should be
    lifted somehow in future...

`parallel_map!()` is used to apply elementwise functions to arbitrary numbers of arguments
using shared-memory parallelism. We do this rather than writing the loops out explicitly
so that `newton_solve!()` and `linear_solve!()` can work for arrays with any combination
of dimensions.

Useful references:
[1] V.A. Mousseau and D.A. Knoll, "Fully Implicit Kinetic Solution of Collisional Plasmas", Journal of Computational Physics 136, 308–323 (1997), https://doi.org/10.1006/jcph.1997.5736.
[2] V.A. Mousseau, "Fully Implicit Kinetic Modelling of Collisional Plasmas", PhD thesis, Idaho National Engineering Laboratory (1996), https://inis.iaea.org/collection/NCLCollectionStore/_Public/27/067/27067141.pdf.
[3] https://en.wikipedia.org/wiki/Generalized_minimal_residual_method
[4] https://www.rikvoorhaar.com/blog/gmres
[5] E. Carson , J. Liesen, Z. Strakoš, "Towards understanding CG and GMRES through examples", Linear Algebra and its Applications 692, 241–291 (2024), https://doi.org/10.1016/j.laa.2024.04.003. 
[6] Q. Zou, "GMRES algorithms over 35 years", Applied Mathematics and Computation 445, 127869 (2023), https://doi.org/10.1016/j.amc.2023.127869
"""
module nonlinear_solvers

export setup_nonlinear_solve, gather_nonlinear_solver_counters!,
       reset_nonlinear_per_stage_counters!, newton_solve!

using ..array_allocation: allocate_float, allocate_shared_float
using ..communication
using ..communication: _anyzv_subblock_synchronize
using ..coordinates: coordinate
using ..debugging: @debug_shared_array
using ..input_structs
using ..jacobian_matrices: create_jacobian_info
using ..looping
using ..timer_utils
using ..type_definitions: mk_float, mk_int

using LinearAlgebra
using MPI
using SparseArrays
using StatsBase: mean

struct nl_solver_info{TH,TV,Tcsg,Tlig,Tprecon,Tpretype,Towned}
    rtol::mk_float
    atol::mk_float
    nonlinear_max_iterations::mk_int
    linear_rtol::mk_float
    linear_atol::mk_float
    linear_restart::mk_int
    linear_max_restarts::mk_int
    H::TH
    c::Tcsg
    s::Tcsg
    g::Tcsg
    V::TV
    linear_initial_guess::Tlig
    n_solves::Base.RefValue{mk_int}
    nonlinear_iterations::Base.RefValue{mk_int}
    linear_iterations::Base.RefValue{mk_int}
    precon_iterations::Base.RefValue{mk_int}
    global_n_solves::Base.RefValue{mk_int}
    global_nonlinear_iterations::Base.RefValue{mk_int}
    global_linear_iterations::Base.RefValue{mk_int}
    global_precon_iterations::Base.RefValue{mk_int}
    solves_since_precon_update::Base.RefValue{mk_int}
    precon_dt::Base.RefValue{mk_float}
    precon_lowerz_vcut_inds::Vector{mk_int}
    precon_upperz_vcut_inds::Vector{mk_int}
    max_nonlinear_iterations_this_step::Base.RefValue{mk_int}
    max_linear_iterations_this_step::Base.RefValue{mk_int}
    total_its_soft_limit::mk_int
    preconditioner_type::Tpretype
    preconditioner_update_interval::mk_int
    preconditioners::Tprecon
    # Ranks in `shared_comm` divide the local entries of each unknown (`local_ranges`).
    # `distributed_comm` is `MPI.COMM_NULL` except on rank 0 of `shared_comm`.
    # `owned_indices` is one tuple entry per unknown array. Each entry is the
    # pair `(local_ranges, owned_ranges)`, fastest index first. `local_ranges`
    # is every slot this rank writes, including the duplicated boundary.
    # `owned_ranges` is the subset of `local_ranges` with that duplicated point
    # omitted on every block except the last.
    shared_comm::MPI.Comm
    distributed_comm::MPI.Comm
    owned_indices::Towned
    # One (sum, count) pair per unknown. Reused by every norm call.
    reduction_buffer::Vector{mk_float}
    full_vector_norm::Bool
end

const local_ranges_index = 1 # index for local points in owned_indices
const owned_ranges_index = 2 # index for points without duplicates in owned_indices

"""

`coords` is a NamedTuple of coordinates corresponding to the dimensions of the variable
that will be solved. Callers pass `coords` slowest-varying first (the reverse of the
array memory layout, which is fastest-varying first).

The nonlinear solver will be called inside a loop over `outer_coords`, so we might need
for example a preconditioner object for each point in that outer loop.
"""
function setup_nonlinear_solve(nl_solver_input, coords, outer_coords=();
                               section_name="nonlinear_solver", serial_solve=false,
                               anysv_region=false, anyzv_region=false,
                               electron_p_pdf_solve=false, preconditioner_type=Val(:none),
                               preconditioners=nothing,
                               full_vector_norm=get(nl_solver_input, :full_vector_norm, false))
    shared_comm, distributed_comm =
        _solver_comms(coords, serial_solve, anysv_region, anyzv_region)

    coord_sizes = Tuple(isa(c, coordinate) ? c.n : c for c ∈ coords)
    outer_coord_sizes = Tuple(isa(c, coordinate) ? c.n : c for c ∈ outer_coords)

    linear_restart = nl_solver_input.linear_restart

    owned_indices = _solver_owned_indices(coords, shared_comm;
                                          electron_p_pdf_solve=electron_p_pdf_solve)
    reduction_buffer = zeros(mk_float, 2 * length(owned_indices))

    n_vcut_inds = 0
    # Rank local variables; only V is a shared allocate when distributed.
    H = allocate_float(linear_restart + 1, linear_restart)
    c = allocate_float(linear_restart + 1)
    s = allocate_float(linear_restart + 1)
    g = allocate_float(linear_restart + 1)
    H .= 0.0
    c .= 0.0
    s .= 0.0
    g .= 0.0
    if shared_comm == MPI.COMM_SELF
        V = allocate_float(reverse(coord_sizes)..., linear_restart+1)
    elseif electron_p_pdf_solve
        V_ppar = allocate_shared_float(coords.z, :linear_restart_plus_one=>linear_restart+1; comm=shared_comm)
        V_pdf = allocate_shared_float(reverse(coords)...,
                                      :linear_restart_plus_one=>linear_restart+1,
                                      comm=shared_comm)
        V = (V_ppar, V_pdf)
        n_vcut_inds = prod(outer_coord_sizes)
    else
        V = allocate_shared_float(reverse(coords)...,
                                  :linear_restart_plus_one=>linear_restart + 1;
                                  comm=shared_comm)
    end
    for i ∈ 1:linear_restart+1
        parallel_map!(() -> 0.0, select_from_V(V, i), owned_indices)
    end
    if shared_comm != MPI.COMM_NULL && shared_comm != MPI.COMM_SELF
        MPI.Barrier(shared_comm)
    end

    linear_initial_guess = zeros(linear_restart)

    return nl_solver_info(mk_float(nl_solver_input.rtol), mk_float(nl_solver_input.atol),
                          nl_solver_input.nonlinear_max_iterations,
                          mk_float(nl_solver_input.linear_rtol),
                          mk_float(nl_solver_input.linear_atol), linear_restart,
                          nl_solver_input.linear_max_restarts, H, c, s, g, V,
                          linear_initial_guess, Ref(0), Ref(0), Ref(0), Ref(0), Ref(0),
                          Ref(0), Ref(0), Ref(0),
                          Ref(nl_solver_input.preconditioner_update_interval),
                          Ref(mk_float(0.0)), zeros(mk_int, n_vcut_inds),
                          zeros(mk_int, n_vcut_inds), Ref(0), Ref(0),
                          nl_solver_input.total_its_soft_limit, preconditioner_type,
                          nl_solver_input.preconditioner_update_interval, preconditioners,
                          shared_comm, distributed_comm, owned_indices,
                          reduction_buffer, full_vector_norm)
end

# `shared_comm` is the communicator whose ranks share each unknown array.
# `distributed_comm` connects the rank-0 leads of `shared_comm` across distributed
# blocks (`MPI.COMM_NULL` on non-zero ranks of `shared_comm` and when not distributed).
function _solver_comms(coords, serial_solve::Bool, anysv_region::Bool, anyzv_region::Bool)
    if serial_solve
        return MPI.COMM_SELF, MPI.COMM_NULL
    elseif anysv_region
        shared_comm = comm_anysv_subblock[]
        return shared_comm, _distributed_comm_for_lead_rank(shared_comm, MPI.COMM_NULL)
    elseif anyzv_region
        shared_comm = comm_anyzv_subblock[]
        return shared_comm, _distributed_comm_for_lead_rank(shared_comm, coords.z.comm)
    else
        shared_comm = comm_block[]
        return shared_comm, _distributed_comm_for_lead_rank(shared_comm, comm_inter_block[])
    end
end

function _distributed_comm_for_lead_rank(shared_comm::MPI.Comm, lead_comm::MPI.Comm)
    if shared_comm == MPI.COMM_NULL || MPI.Comm_rank(shared_comm) != 0
        return MPI.COMM_NULL
    end
    return lead_comm
end

"""
    _solver_owned_indices(coords, shared_comm; electron_p_pdf_solve)

Flat index ranges for this rank. Callers pass `coords` slowest-varying first. The
ranges index array memory, which is fastest-varying first. The entries are divided across
`shared_comm`. A rank with `MPI.COMM_NULL`, or a rank in the communicator that receives
no entries, gets an empty list. Membership of `shared_comm` is unchanged.

For each unknown the pair `(local_ranges, owned_ranges)` is returned. The full local
array (`1:prod(dim_sizes)`) is partitioned once across `shared_comm` to form
`local_ranges`, which is every slot this rank writes (including the duplicated upper
block-boundary point). `owned_ranges` is formed by intersecting `local_ranges` with the
non-duplicated block points (omitting the upper boundary point on every block except the
last along each distributed coordinate), guaranteeing `owned_ranges ⊆ local_ranges` on
every rank without shifting the inter-rank cut. Consecutive indices are stored as one
range.

A single unknown returns a one-element tuple whose element is that pair. Implicit `vpa`
advection, Braginskii conduction, Fokker-Planck, and the implicit ion advance use this
form.

`electron_p_pdf_solve` returns two pairs, pressure then the distribution. Pressure uses
`coords.z` only. The distribution uses every dimension in `coords`. The pairs are split
separately, so each unknown keeps its own norm.
"""
function _solver_owned_indices(coords, shared_comm::MPI.Comm; electron_p_pdf_solve::Bool)
    if !electron_p_pdf_solve
        dim_sizes, exclude_last = _sizes_and_excluded(coords)
        return (_local_and_owned_ranges(dim_sizes, exclude_last, shared_comm),)
    end
    z = coords.z
    pressure_ranges = _local_and_owned_ranges(
        (Int(z.n),), (_dimension_excludes_boundary(z),), shared_comm)
    pdf_sizes, pdf_exclude = _sizes_and_excluded(coords)
    pdf_ranges = _local_and_owned_ranges(pdf_sizes, pdf_exclude, shared_comm)
    return pressure_ranges, pdf_ranges
end

# Split the full local array first, then intersect with `owned_block` so every rank's
# `owned_ranges` is a subset of its `local_ranges`.
function _local_and_owned_ranges(dim_sizes::NTuple{N,Int}, exclude_last::NTuple{N,Bool},
                                 shared_comm::MPI.Comm) where {N}
    local_block = _block_owned_ranges(dim_sizes, ntuple(_ -> false, N))
    local_ranges = _split_ranges_across_comm(local_block, shared_comm)
    owned_block = _block_owned_ranges(dim_sizes, exclude_last)
    owned_ranges = _intersect_ranges(local_ranges, owned_block)
    return (local_ranges, owned_ranges) # index with local_ranges_index, owned_ranges_index
end

# `coords` is ordered slowest-varying first. Array memory is fastest-varying first.
function _sizes_and_excluded(coords::NamedTuple{Names}) where {Names}
    N = length(Names)
    dim_sizes = ntuple(Val(N)) do i
        c = coords[i]
        Int(isa(c, coordinate) ? c.n : c)
    end
    exclude_last = ntuple(Val(N)) do i
        _dimension_excludes_boundary(coords[i])
    end
    return reverse(dim_sizes), reverse(exclude_last)
end

function _dimension_excludes_boundary(c)
    return c isa coordinate && c.nrank > 1 && c.irank < c.nrank - 1
end

# Contiguous flat-index ranges for the points owned by this distributed-memory
# block. The last local index of a distributed dimension is the copy of the
# next block's boundary point, so it is omitted.
function _block_owned_ranges(dim_sizes::NTuple{N,Int}, exclude_last::NTuple{N,Bool}) where {N}
    owned_sizes = ntuple(d -> exclude_last[d] ? dim_sizes[d] - 1 : dim_sizes[d], N)
    if any(s -> s <= 0, owned_sizes)
        return UnitRange{Int}[]
    end
    if N == 1
        return UnitRange{Int}[1:owned_sizes[1]]
    end
    inner = owned_sizes[1]
    ranges = UnitRange{Int}[]
    outer_index = ones(Int, N - 1)
    nouter = prod(owned_sizes[2:end])
    sizehint!(ranges, nouter)
    for _ in 1:nouter
        offset = 0
        stride = dim_sizes[1]
        for j in 1:(N - 1)
            offset += (outer_index[j] - 1) * stride
            stride *= dim_sizes[j + 1]
        end
        push!(ranges, (offset + 1):(offset + inner))
        for j in 1:(N - 1)
            outer_index[j] += 1
            if outer_index[j] <= owned_sizes[j + 1]
                break
            end
            outer_index[j] = 1
        end
    end
    return _coalesce_ranges(ranges)
end

function _coalesce_ranges(ranges::Vector{UnitRange{Int}})
    isempty(ranges) && return ranges
    coalesced = UnitRange{Int}[]
    start_ind = first(first(ranges))
    stop_ind = last(first(ranges))
    for r in @view ranges[2:end]
        if first(r) == stop_ind + 1
            stop_ind = last(r)
        else
            push!(coalesced, start_ind:stop_ind)
            start_ind = first(r)
            stop_ind = last(r)
        end
    end
    push!(coalesced, start_ind:stop_ind)
    return coalesced
end

function _split_ranges_across_comm(block_ranges::Vector{UnitRange{Int}}, comm::MPI.Comm)
    if comm == MPI.COMM_NULL
        return UnitRange{Int}[]
    end
    n_owned = sum(length, block_ranges; init=0)
    mine = _chunk_of_owned_count(n_owned, MPI.Comm_rank(comm), MPI.Comm_size(comm))
    return _slice_owned_ranges(block_ranges, mine)
end

# Cut the same as `get_local_range`: extra points go to the higher ranks, so this
# rank reads the entries it writes.
function _chunk_of_owned_count(n_owned::Int, rank::Int, nrank::Int)
    base, remainder = divrem(n_owned, nrank)
    n_low = nrank - remainder
    if rank < n_low
        start_ind = rank * base + 1
        return start_ind:(start_ind + base - 1)
    else
        start_ind = n_low * base + (rank - n_low) * (base + 1) + 1
        return start_ind:(start_ind + base)
    end
end

function _slice_owned_ranges(block_ranges::Vector{UnitRange{Int}}, mine::UnitRange{Int})
    if isempty(mine)
        return UnitRange{Int}[]
    end
    sliced = UnitRange{Int}[]
    position = 1
    for r in block_ranges
        rlen = length(r)
        range_end = position + rlen - 1
        overlap_start = max(first(mine), position)
        overlap_stop = min(last(mine), range_end)
        if overlap_start <= overlap_stop
            shift = overlap_start - position
            len = overlap_stop - overlap_start + 1
            start_ind = first(r) + shift
            push!(sliced, start_ind:(start_ind + len - 1))
        end
        position = range_end + 1
        if position > last(mine)
            break
        end
    end
    return sliced
end

function _intersect_ranges(a::Vector{UnitRange{Int}}, b::Vector{UnitRange{Int}})
    result = UnitRange{Int}[]
    ia = 1
    ib = 1
    while ia <= length(a) && ib <= length(b)
        ra = a[ia]
        rb = b[ib]
        start_ind = max(first(ra), first(rb))
        stop_ind = min(last(ra), last(rb))
        if start_ind <= stop_ind
            push!(result, start_ind:stop_ind)
        end
        if last(ra) < last(rb)
            ia += 1
        else
            ib += 1
        end
    end
    return result
end

"""
    reset_nonlinear_per_stage_counters!(nl_solver_params::Union{nl_solver_info,Nothing})

Reset the counters that hold per-step totals or maximums in `nl_solver_params`.

Also increment `nl_solver_params.stage_counter[]`.
"""
function reset_nonlinear_per_stage_counters!(nl_solver_params::Union{nl_solver_info,Nothing})
    if nl_solver_params === nothing
        return nothing
    end

    nl_solver_params.max_nonlinear_iterations_this_step[] = 0
    nl_solver_params.max_linear_iterations_this_step[] = 0

    # Also increment the stage counter
    nl_solver_params.solves_since_precon_update[] += 1

    return nothing
end

"""
    gather_nonlinear_solver_counters!(nl_solver_params)

Where necessary, gather the iteration counters for the nonlinear solvers.

Where each solve runs in parallel using all processes, this is unnecessary as the count on
each process already represents the global count. Where each solve uses only a subset of
processes, the counters from different solves need to be added together to get the global
total.
"""
@timeit_debug global_timer gather_nonlinear_solver_counters!(nl_solver_params) = begin
    if nl_solver_params.ion_advance !== nothing
        # Solve runs in parallel on all processes, so no need to collect here
        nl_solver_params.ion_advance.global_n_solves[] = nl_solver_params.ion_advance.n_solves[]
        nl_solver_params.ion_advance.global_nonlinear_iterations[] = nl_solver_params.ion_advance.nonlinear_iterations[]
        nl_solver_params.ion_advance.global_linear_iterations[] = nl_solver_params.ion_advance.linear_iterations[]
        nl_solver_params.ion_advance.global_precon_iterations[] = nl_solver_params.ion_advance.precon_iterations[]
    end
    if nl_solver_params.vpa_advection !== nothing
        # Solves are run in serial on separate processes, so need a global Allreduce
        @timeit_debug global_timer "MPI.Allreduce! comm_world" MPI.Allreduce!(nl_solver_params.vpa_advection.n_solves[], +, comm_world)
        @timeit_debug global_timer "MPI.Allreduce! comm_world" MPI.Allreduce!(nl_solver_params.vpa_advection.nonlinear_iterations[], +, comm_world)
        @timeit_debug global_timer "MPI.Allreduce! comm_world" MPI.Allreduce!(nl_solver_params.vpa_advection.linear_iterations[], +, comm_world)
        @timeit_debug global_timer "MPI.Allreduce! comm_world" MPI.Allreduce!(nl_solver_params.vpa_advection.precon_iterations[], +, comm_world)
    end
end

"""
    newton_solve!(x, rhs_func!, residual, delta_x, rhs_delta, w, nl_solver_params;
                  left_preconditioner=identity, right_preconditioner=identity, coords)

`x` is the initial guess at the solution, and is overwritten by the result of the Newton
solve.

`rhs_func!(residual, x)` is the function we are trying to find a solution of. It calculates
```math
\\mathtt{residual} = F(\\mathtt{x})
```
where we are trying to solve \$F(x)=0\$.

`residual`, `delta_x`, `rhs_delta` and `w` are buffer arrays, with the same size as `x`,
used internally.

`left_preconditioner` or `right_preconditioner` apply preconditioning. They should be
passed a function that solves \$P.x = b\$ where \$P\$ is the preconditioner matrix, \$b\$
is given by the values passed to the function as the argument, and the result \$x\$ is
returned by overwriting the argument.

`coords` is a NamedTuple containing the `coordinate` structs corresponding to each
dimension in `x`.


Tolerances
----------

Note that the meaning of the relative tolerance `rtol` and absolute tolerance `atol` is
very different for the outer Newton iteration and the inner GMRES iteration.

For the outer Newton iteration the residual \$R(x^n)\$ measures the departure of the
system from the solution (at each grid point). Its size can be compared to the size of the
solution `x`, so it makes sense to define an `error norm' for \$R(x^n)\$ as
```math
E(x^n) = \\left\\lVert \\frac{R(x^n)}{\\mathtt{rtol} x^n \\mathtt{atol}} \\right\\rVert_2
```
where \$\\left\\lVert \\cdot \\right\\rVert\$ is the 'L2 norm' (square-root of sum of
squares). We can further try to define a grid-size independent error norm by dividing out
the number of grid points to get a root-mean-square (RMS) error rather than an L2 norm.
```math
E_{\\mathrm{RMS}}(x^n) = \\sqrt{ \\frac{1}{N} \\sum_i \\frac{R(x^n)_i}{\\mathtt{rtol} x^n_i \\mathtt{atol}} }
```
where \$N\$ is the total number of grid points.

In contrast, GMRES is constructed to minimise the L2 norm of \$r_k = b - A\\cdot x_k\$
where GMRES is solving the linear system \$A\\cdot x = b\$, \$x_k\$ is the approximation
to the solution \$x\$ at the \$k\$'th iteration and \$r_k\$ is the residual at the
\$k\$'th iteration. There is no flexibility to measure error relative to \$x\$ in any
sense. For GMRES, a `relative tolerance' is relative to the residual of the
right-hand-side \$b\$, which is the first iterate \$x_0\$ (when no initial guess is
given). [Where a non-zero initial guess is given it might be better to use a different
stopping criterion, see Carson et al. section 3.8.]. The stopping criterion for the GMRES
iteration is therefore
```
\\left\\lVert r_k \\right\\rVert < \\max(\\mathtt{linear\\_rtol} \\left\\lVert r_0 \\right\\rVert, \\mathtt{linear\\_atol}) = \\max(\\mathtt{linear\\_rtol} \\left\\lVert b \\right\\rVert, \\mathtt{linear\\_atol})
```
As the GMRES solve is only used to get the right `direction' for the next Newton step, it
is not necessary to have a very tight `linear_rtol` for the GMRES solve.
"""
@timeit global_timer newton_solve!(
                         x, residual_func!::F, residual, delta_x, rhs_delta, v, w,
                         nl_solver_params; left_preconditioner::LP=identity,
                         right_preconditioner::RP=identity,
                         recalculate_preconditioner::RECALC=nothing,
                         coords=nothing) where {F,LP,RP,RECALC} = begin
    rtol = nl_solver_params.rtol
    atol = nl_solver_params.atol
    owned_indices = nl_solver_params.owned_indices
    shared_comm = nl_solver_params.shared_comm

    norm_params = (nl_solver_params.owned_indices,
                   nl_solver_params.shared_comm,
                   nl_solver_params.distributed_comm,
                   nl_solver_params.reduction_buffer, rtol, atol, x,
                   nl_solver_params.full_vector_norm)

    residual_func!(residual, x)
    # The residual writes with @loop_*. Synchronize before the owned-range read.
    _sync_shared(shared_comm)
    residual_norm = distributed_norm(residual, norm_params...)
    counter = 0
    linear_counter = 0

    # Would need this if delta_x was not set to zero within the Newton iteration loop
    # below.
    #parallel_map!(()->0.0, delta_x, owned_indices)

    close_counter = -1
    close_linear_counter = -1
    success = true
    previous_residual_norm = residual_norm
    old_precon_iterations = nl_solver_params.precon_iterations[]
    while (counter < 1 && residual_norm > 1.0e-8) || residual_norm > 1.0
        counter += 1
        #println("\nNewton ", counter)

        # Solve (approximately?):
        #   J δx = -RHS(x)
        parallel_map!(()->0.0, delta_x, owned_indices)
        linear_its = linear_solve!(x, residual_func!, residual, delta_x, v, w,
                                   norm_params;
                                   rtol=nl_solver_params.linear_rtol,
                                   atol=nl_solver_params.linear_atol,
                                   restart=nl_solver_params.linear_restart,
                                   max_restarts=nl_solver_params.linear_max_restarts,
                                   left_preconditioner=left_preconditioner,
                                   right_preconditioner=right_preconditioner,
                                   H=nl_solver_params.H, c=nl_solver_params.c,
                                   s=nl_solver_params.s, g=nl_solver_params.g,
                                   V=nl_solver_params.V, rhs_delta=rhs_delta,
                                   initial_guess=nl_solver_params.linear_initial_guess,
                                   initial_delta_x_is_zero=true)
        linear_counter += linear_its
        # delta_x was written on local_ranges (or by the preconditioner). Synchronize.
        _sync_shared(shared_comm)

        # If the residual does not decrease, we will do a line search to find an update
        # that does decrease the residual. The value of `x` is used to define the
        # normalisation value with rtol that is used to calculate the residual, so do not
        # want to update it until the line search is completed (otherwise the norm changes
        # during the line search, which might make it fail to converge). So calculate the
        # updated value in the buffer `w` until the line search is completed, and only
        # then copy it into `x`.
        parallel_map!((x) -> x, w, x, owned_indices)
        parallel_map!((x,delta_x) -> x + delta_x, w, x, delta_x, owned_indices)
        # w was written on local_ranges. The residual reads it with @loop_*.
        _sync_shared(shared_comm)
        residual_func!(residual, w)
        # The residual writes with @loop_*. Synchronize before continuing.
        _sync_shared(shared_comm)

        # For the Newton iteration, we want the norm divided by the (sqrt of the) number
        # of grid points, so we can use a tolerance that is independent of the size of the
        # grid. This is unlike the norms needed in `linear_solve!()`.
        residual_norm = distributed_norm(residual, norm_params...)
        if isnan(residual_norm)
            error("NaN in Newton iteration at iteration $counter")
        end
#        if residual_norm > previous_residual_norm
#            # Do a line search between x and x+delta_x to try to find an update that does
#            # decrease residual_norm
#            s = 0.5
#            while s > 1.0e-2
#                parallel_map!(solver_type, (x,delta_x,s) -> x + s * delta_x, w, x, delta_x, s)
#                residual_func!(residual, x)
#                residual_norm = distributed_norm(solver_type, residual, norm_params...)
#                if residual_norm ≤ previous_residual_norm
#                    break
#                end
#                s *= 0.5
#            end
#            println("line search s ", s)
#
#            #if residual_norm > previous_residual_norm
#            #    # Failed to find a point that decreases the residual, so try a negative
#            #    # step
#            #    s = -1.0e-5
#            #    parallel_map!(solver_type, (x,delta_x,s) -> x + s * delta_x, w, x, delta_x, s)
#            #    residual_func!(residual, x)
#            #    residual_norm = distributed_norm(solver_type, residual, norm_params...)
#            #    if residual_norm > previous_residual_norm
#            #        # That didn't work either, so just take the full step and hope for
#            #        # convergence later
#            #        parallel_map!(solver_type, (x,delta_x,s) -> x + s * delta_x, w, x, delta_x, s)
#            #        residual_func!(residual, x)
#            #        residual_norm = distributed_norm(solver_type, residual, norm_params...)
#            #    end
#            #end
#            if residual_norm > previous_residual_norm
#                # Line search didn't work, so just take the full step and hope for
#                # convergence later
#                parallel_map!(solver_type, (x,delta_x,s) -> x + s * delta_x, w, x, delta_x, s)
#                residual_func!(residual, x)
#                residual_norm = distributed_norm(solver_type, residual, norm_params...)
#            end
#        end
        parallel_map!((w) -> w, x, w, owned_indices)
        # x was written on local_ranges. Synchronize for preconditioner updates and the caller.
        _sync_shared(shared_comm)
        previous_residual_norm = residual_norm

        if recalculate_preconditioner !== nothing && counter % nl_solver_params.preconditioner_update_interval == 0
            # Have taken a large number of Newton iterations already - convergence must be
            # slow, so try updating the preconditioner.
            recalculate_preconditioner()
        end

        #println("Newton residual ", residual_norm, " ", linear_its, " $rtol $atol")

        if residual_norm < 0.1/rtol && close_counter < 0 && close_linear_counter < 0
            close_counter = counter
            close_linear_counter = linear_counter
        end

        if counter > nl_solver_params.nonlinear_max_iterations
            println("maximum iteration limit reached")
            success = false
            break
        end
    end
    if counter == 0
        _sync_shared(shared_comm)
    end
    nl_solver_params.n_solves[] += 1
    nl_solver_params.nonlinear_iterations[] += counter
    nl_solver_params.linear_iterations[] += linear_counter
    nl_solver_params.max_nonlinear_iterations_this_step[] =
        max(counter, nl_solver_params.max_nonlinear_iterations_this_step[])
    nl_solver_params.max_linear_iterations_this_step[] =
        max(linear_counter, nl_solver_params.max_linear_iterations_this_step[])
#    println("Newton iterations: ", counter)
#    println("Final residual: ", residual_norm)
#    println("Total linear iterations: ", linear_counter)
#    println("Linear iterations per Newton: ", linear_counter / counter)
#    precon_count = nl_solver_params.precon_iterations[] - old_precon_iterations
#    println("Total precon iterations: ", precon_count)
#    println("Precon iterations per linear: ", precon_count / linear_counter)
#
#    println("Newton iterations after close: ", counter - close_counter)
#    println("Total linear iterations after close: ", linear_counter - close_linear_counter)
#    println("Linear iterations per Newton after close: ", (linear_counter - close_linear_counter) / (counter - close_counter))
#    println()

    return success
end

# Weighted RMS of the scaled entries. `residual` is the vector being measured
# and `x` is the current Newton state, used as the scale `rtol*abs(x)+atol`.
# Works with a tuple of arrays; a bare array is wrapped as a one-element list.
# When `full_vector_norm == false` (default), each unknown is divided by its own
# point count, then those mean squares are averaged so that multiple unknowns
# have equal weight. When `full_vector_norm == true`, all owned points across all
# unknowns are summed together and divided by the total owned point count.
@timeit_debug global_timer distributed_norm(residual, owned_indices,
                          shared_comm, distributed_comm, buffer, rtol, atol, x,
                          full_vector_norm::Bool=false) = begin
    mean = _weighted_mean(residual, residual, owned_indices, shared_comm,
                          distributed_comm, buffer, rtol, atol, x,
                          full_vector_norm) do rk, _, i
        rk[i] * rk[i]
    end
    return sqrt(mean)
end

# Same mean as `distributed_norm`, with the product of `v` and `w` as the weight.
@timeit_debug global_timer distributed_dot(v, w, owned_indices,
                          shared_comm, distributed_comm, buffer, rtol, atol, x,
                          full_vector_norm::Bool=false) = begin
    return _weighted_mean(v, w, owned_indices, shared_comm,
                          distributed_comm, buffer, rtol, atol, x,
                          full_vector_norm) do vk, wk, i
        vk[i] * wk[i]
    end
end

# Mean of `weight(vk, wk, i) / scale^2`, with `scale = rtol*abs(x)+atol`.
# Slot 2k-1 is the sum of those weights and slot 2k is the point count.
# `_reduce_sums!` runs before the division. When `full_vector_norm` is false,
# the per-unknown means are averaged; when true, total sum is divided by total count.
# `map` over the tuples unrolls via tuple recursion so heterogeneous unknown tuples
# (e.g. 1D pressure + 3D pdf) specialize and SIMD-vectorize per unknown.
@inline function _weighted_mean(weight::F, v, w, owned_indices, shared_comm,
                                distributed_comm, buffer, rtol, atol, x,
                                full_vector_norm::Bool=false) where {F}
    vs = v isa AbstractArray ? (v,) : v
    ws = w isa AbstractArray ? (w,) : w
    xs = x isa AbstractArray ? (x,) : x
    n = length(xs)
    totals = map(owned_indices, vs, ws, xs) do oi, vk, wk, xk
        s = zero(mk_float)
        c = 0
        # Owned ranges omit the boundary point duplicated on the next block.
        for r in oi[owned_ranges_index]
            for i in r
                scale = rtol * abs(xk[i]) + atol
                s += weight(vk, wk, i) / scale^2
            end
            c += length(r)
        end
        return (s, c)
    end
    for k in 1:n
        buffer[2k - 1], buffer[2k] = totals[k]
    end
    _reduce_sums!(buffer, shared_comm, distributed_comm)
    if full_vector_norm
        # Every grid point weighted equally across the whole solution vector.
        total_sum = 0.0
        total_count = 0.0
        for k in 1:n
            total_sum += buffer[2k - 1]
            total_count += buffer[2k]
        end
        return total_sum / total_count
    else
        # Unknowns weighted equally.
        mean = 0.0
        for k in 1:n
            mean += buffer[2k - 1] / buffer[2k]
        end
        return mean / n
    end
end

# Synchronize ranks in `shared_comm` when crossing between `local_ranges` and `@loop_*`
# (or preconditioners) so shared-array writes are visible before the next read or overwrite.
function _sync_shared(shared_comm::MPI.Comm)
    if shared_comm == MPI.COMM_NULL || shared_comm == MPI.COMM_SELF
        return nothing
    end
    @debug_shared_array begin
        if shared_comm == comm_anyzv_subblock[]
            @_anyzv_subblock_synchronize()
            return nothing
        elseif shared_comm == comm_anysv_subblock[]
            @_anysv_subblock_synchronize()
            return nothing
        elseif shared_comm == comm_block[]
            @_block_synchronize()
            return nothing
        end
    end
    MPI.Barrier(shared_comm)
    return nothing
end

# Sum-reduce the (sum, count) pairs to every rank.
function _reduce_sums!(buffer, shared_comm, distributed_comm)
    if shared_comm != MPI.COMM_NULL
        # The total is present on rank 0 only until the broadcast below.
        MPI.Reduce!(buffer, +, shared_comm)
    end
    if distributed_comm != MPI.COMM_NULL
        # Sum across blocks. Allreduce leaves the total on every rank in this comm.
        MPI.Allreduce!(buffer, +, distributed_comm)
    end
    if shared_comm != MPI.COMM_NULL
        # Copy rank 0's buffer to the other ranks in the shared group.
        MPI.Bcast!(buffer, shared_comm; root=0)
    end
    return buffer
end

# 0 input parallel_map!
@inline function parallel_map!(func::F, result, owned_indices) where {F}
    value = func()
    results = result isa AbstractArray ? (result,) : result
    map(owned_indices, results) do oi, field
        for r in oi[local_ranges_index]
            for idx in r
                field[idx] = value
            end
        end
        return nothing
    end
    return nothing
end

# 1 input parallel_map!
@inline function parallel_map!(func::F, result, x1, owned_indices) where {F}
    results = result isa AbstractArray ? (result,) : result
    x1s = x1 isa AbstractArray ? (x1,) : x1
    map(owned_indices, results, x1s) do oi, field, in1
        for r in oi[local_ranges_index]
            for idx in r
                field[idx] = func(in1[idx])
            end
        end
        return nothing
    end
    return nothing
end

# 2 input parallel_map!
@inline function parallel_map!(func::F, result, x1, x2, owned_indices) where {F}
    results = result isa AbstractArray ? (result,) : result
    x1s = x1 isa AbstractArray ? (x1,) : x1
    if x2 isa AbstractArray || x2 isa Tuple
        x2s = x2 isa AbstractArray ? (x2,) : x2
        map(owned_indices, results, x1s, x2s) do oi, field, in1, in2
            for r in oi[local_ranges_index]
                for idx in r
                    field[idx] = func(in1[idx], in2[idx])
                end
            end
            return nothing
        end
    else
        map(owned_indices, results, x1s) do oi, field, in1
            for r in oi[local_ranges_index]
                for idx in r
                    field[idx] = func(in1[idx], x2)
                end
            end
            return nothing
        end
    end
    return nothing
end


@timeit_debug global_timer parallel_delta_x_calc!(delta_x, V, y, owned_indices) = begin
    for iy ∈ 1:length(y)
        y_iy = y[iy]
        parallel_map!((dx, v) -> dx + y_iy * v, delta_x, delta_x, select_from_V(V, iy),
                     owned_indices)
    end
    return nothing
end

# Utility function for neatness handling that V may be an array or a Tuple of arrays
function select_from_V(V::Tuple, i)
    return map(Vpart -> selectdim(Vpart, ndims(Vpart), i), V)
end
function select_from_V(V, i)
    return selectdim(V,ndims(V),i)
end

@kwdef struct nonlinear_solvers_approximate_Jacobian_vector_product!{F,LP,RP,Tx,Tr0,Tdel,Tidx}
    residual_func!::F
    left_preconditioner::LP
    right_preconditioner::RP
    x::Tx
    residual0::Tr0
    rhs_delta::Tdel
    Jv_scale_factor::mk_float
    inv_Jv_scale_factor::mk_float
    owned_indices::Tidx
    shared_comm::MPI.Comm
end

function (approx_Jv::nonlinear_solvers_approximate_Jacobian_vector_product!)(
             v, skip_first_precon::Bool=false)
    if !skip_first_precon && approx_Jv.right_preconditioner !== identity
        _sync_shared(approx_Jv.shared_comm)
        approx_Jv.right_preconditioner(v)
        _sync_shared(approx_Jv.shared_comm)
    end

    owned_indices = approx_Jv.owned_indices
    rhs_delta = approx_Jv.rhs_delta

    parallel_map!((x,v) -> x + approx_Jv.Jv_scale_factor * v, v, approx_Jv.x, v,
                 owned_indices)
    # The probe was written on local_ranges. The residual reads it with @loop_*.
    _sync_shared(approx_Jv.shared_comm)
    approx_Jv.residual_func!(rhs_delta, v; krylov=true)
    # The residual wrote with @loop_*. Synchronize before the local_ranges read.
    _sync_shared(approx_Jv.shared_comm)
    parallel_map!((rhs_delta, residual0) -> (rhs_delta - residual0) * approx_Jv.inv_Jv_scale_factor,
                 v, rhs_delta, approx_Jv.residual0, owned_indices)
    if approx_Jv.left_preconditioner !== identity
        _sync_shared(approx_Jv.shared_comm)
        approx_Jv.left_preconditioner(v)
        _sync_shared(approx_Jv.shared_comm)
    end
    return v
end

"""
Apply the GMRES algorithm to solve the 'linear problem' J.δx^n = R(x^n), which is needed
at each step of the outer Newton iteration (in `newton_solve!()`).

Uses Givens rotations to reduce the upper Hessenberg matrix to an upper triangular form,
which allows conveniently finding the residual at each step, and computing the final
solution, without calculating a least-squares minimisation at each step. See 'algorithm 2
MGS-GMRES' in Zou (2023) [https://doi.org/10.1016/j.amc.2023.127869].
"""
@timeit global_timer linear_solve!(
                         x, residual_func!::F, residual0, delta_x, v, w,
                         norm_params; rtol, atol, restart, max_restarts,
                         left_preconditioner::LP, right_preconditioner::RP, H, c, s, g, V,
                         rhs_delta, initial_guess, initial_delta_x_is_zero) where {F,LP,RP} = begin
    # Solve (approximately?):
    #   J δx = residual0

    Jv_scale_factor = 1.0e3
    inv_Jv_scale_factor = 1.0 / Jv_scale_factor

    # The vectors `v` that are passed to this function will be normalised so that
    # `distributed_norm(v) == 1.0`. `distributed_norm()` is defined - including the
    # relative and absolute tolerances from the Newton iteration - so that a vector with a
    # norm of 1.0 is 'small' in the sense that a vector with a norm of 1.0 is small enough
    # relative to `x` to consider the iteration converged. This means that `x+v` would be
    # very close to `x`, so R(x+v)-R(x) would be likely to be badly affected by rounding
    # errors, because `v` is so small, relative to `x`. We actually want to multiply `v`
    # by a large number `Jv_scale_factor` (in constrast to the small `epsilon` in the
    # 'usual' case where the norm does not include either reative or absolute tolerance)
    # to ensure that we get a reasonable estimate of J.v.
    owned_indices = norm_params[1]
    shared_comm = norm_params[2]
    approximate_Jacobian_vector_product! =
        nonlinear_solvers_approximate_Jacobian_vector_product!(;
            residual_func!, left_preconditioner, right_preconditioner, x, residual0,
            rhs_delta, Jv_scale_factor, inv_Jv_scale_factor, owned_indices,
            shared_comm)

    parallel_map!((delta_x) -> delta_x, v, delta_x, owned_indices)

    # This function transforms the data stored in 'v' from δx to ≈J.δx
    # If initial δx is all-zero, we can skip a right-preconditioner evaluation because it
    # would just transform all-zero to all-zero.
    approximate_Jacobian_vector_product!(v, initial_delta_x_is_zero)

    # To start with we use 'w' as a buffer to make a copy of residual0 to which we can
    # apply the left-preconditioner.
    parallel_map!((residual_0) -> residual_0, w, residual0, owned_indices)
    if left_preconditioner !== identity
        _sync_shared(shared_comm)
        left_preconditioner(w)
        _sync_shared(shared_comm)
    end

    # Now we actually set 'w' as the first Krylov vector, and normalise it.
    parallel_map!((residual0, v) -> -residual0 - v, w, w, v, owned_indices)
    beta = distributed_norm(w, norm_params...)
    parallel_map!((w,beta) -> w/beta, select_from_V(V, 1), w, beta, owned_indices)
    g[1] = beta

    # Set tolerance for GMRES iteration to rtol times the initial residual, unless this is
    # so small that it is smaller than atol, in which case use atol instead.
    tol = max(rtol * beta, atol)

    lsq_result = nothing
    residual = Inf
    counter = 0
    restart_counter = 1
    while true
        inner_counter = 0
        for i ∈ 1:restart
            inner_counter = i
            counter += 1
            #println("Linear ", counter)

            # Compute next Krylov vector
            parallel_map!((V) -> V, w, select_from_V(V, i), owned_indices)
            approximate_Jacobian_vector_product!(w)

            # Gram-Schmidt orthogonalization
            for j ∈ 1:i
                w_dot_Vj = distributed_dot(w, select_from_V(V, j), norm_params...)
                H[j,i] = w_dot_Vj
                parallel_map!((w, V) -> w - H[j,i] * V, w, w, select_from_V(V, j), owned_indices)
            end
            norm_w = distributed_norm(w, norm_params...)
            H[i+1,i] = norm_w
            parallel_map!((w) -> w / H[i+1,i], select_from_V(V, i+1), w, owned_indices)

            for j ∈ 1:i-1
                gamma = c[j] * H[j,i] + s[j] * H[j+1,i]
                H[j+1,i] = -s[j] * H[j,i] + c[j] * H[j+1,i]
                H[j,i] = gamma
            end
            delta = sqrt(H[i,i]^2 + H[i+1,i]^2)
            s[i] = H[i+1,i] / delta
            c[i] = H[i,i] / delta
            H[i,i] = c[i] * H[i,i] + s[i] * H[i+1,i]
            H[i+1,i] = 0
            g[i+1] = -s[i] * g[i]
            g[i] = c[i] * g[i]
            residual = abs(g[i+1])

            if residual < tol
                break
            end
        end
        i = inner_counter

        # Update initial guess to restart
        #################################

        @views y = H[1:i,1:i] \ g[1:i]

        # The following calculates
        #    delta_x .= delta_x .+ sum(y[i] .* V[:,i] for i ∈ 1:length(y))
        parallel_delta_x_calc!(delta_x, V, y, owned_indices)
        if right_preconditioner !== identity
            _sync_shared(shared_comm)
            right_preconditioner(delta_x)
        end

        if residual < tol || restart_counter > max_restarts
            break
        end

        restart_counter += 1

        # Store J.delta_x in the variable delta_x, to use it to calculate the new first
        # Krylov vector v/beta.
        _sync_shared(shared_comm)
        parallel_map!((delta_x) -> delta_x, v, delta_x, owned_indices)
        approximate_Jacobian_vector_product!(v)

        # Note residual0 has already had the left_preconditioner() applied to it.
        parallel_map!((residual0, v) -> -residual0 - v, v, residual0, v, owned_indices)
        beta = distributed_norm(v, norm_params...)
        for i ∈ 2:length(y)
            parallel_map!(() -> 0.0, select_from_V(V, i), owned_indices)
        end
        parallel_map!((v,beta) -> v/beta, select_from_V(V, 1), v, beta, owned_indices)
    end

    return counter
end

end
