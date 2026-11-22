export coupled_systems

# Made with help from Claude Sonnet 5.5 High

using DynamicalSystems

# All coupling info lives in the parameter container, so it can be changed later
# with `set_parameter!(cds, :k, 0.3)`, `set_parameter!(cds, :τ, [1.0, 2.0])`, etc.
mutable struct CoupledParams{F, P, T <: AbstractVector, A <: AbstractMatrix, J}
    f::F                    # in-place rule of the original system
    p::P                    # parameters of the original system
    τ::T                    # timescale of each copy
    k::Float64              # coupling strength
    adj::A                  # adjacency matrix (N×N), adj[i,m] = 1 → m drives i
    js::J                   # coupled dimensions
    D::Int                  # dimension of the original system
    scale_coupling::Bool    # whether to also multiply the coupling term by τ[i]
end

# These three functions are needed for `set_parameters!` to work (GUI needs that)
Base.length(p::CoupledParams) = nfields(p)
Base.eltype(::Type{<:CoupledParams}) = Pair{Symbol, Any}
function Base.iterate(p::CoupledParams, i::Int = 1)
    i > nfields(p) && return nothing
    T = typeof(p)
    return (fieldname(T, i) => getfield(p, i), i + 1)
end

function coupled_rule!(du, u, p::CoupledParams, t)
    D = p.D
    N = length(p.τ)
    for i in 1:N # loop over coupled systems
        r = ((i - 1) * D + 1):(i * D)
        dui = view(du, r)
        p.f(dui, view(u, r), p.p, t)
        dui .*= p.τ[i]
        c = p.scale_coupling ? p.τ[i] : one(eltype(du))
        for j in p.js
            n = (i - 1) * D + j
            s = zero(eltype(du))
            for m in 1:N
                s += p.adj[i, m] * (u[(m - 1) * D + j] - u[n])
            end
            du[n] += c * p.k * s
        end
    end
    return nothing
end

"""
    couple_systems(ds::CoupledODEs, τ, js; k = 1.0, adjacency = nothing,
                   scale_coupling = false, u0s = nothing, diffeq = ds.diffeq)

Couple `length(τ)` identical copies of the in-place system `ds` into one bigger
in-place `CoupledODEs`. Copy `i` evolves as

    u̇ᵢ = τᵢ f(uᵢ, p, t)  +  k Σₘ Aᵢₘ (uₘ[j] − uᵢ[j])   (for every j ∈ js)

Each copy shares the same internal parameters as `ds`.

## Arguments

- `τ`: vector containing the timescale of each copy (multiplies that copy's RHS).
- `js`: dimension of  (Int) or dimensions (collection) `ds` along which the diffusive coupling acts.
- `k`: coupling strength.
- `adjacency = nothing`: `N×N` matrix deciding which identical copies of `ds` couple to which other.
  If `adj[i, m] = 1` then m drives i. The default is all-to-all without self-coupling.
- `scale_coupling = false`: if `true`, the coupling term is multiplied by `τᵢ` as well,
  i.e. `u̇ᵢ = τᵢ (f + coupling)`.
- `u0s`: vector of initial states, one per copy (default: copies of `current_state(ds)`).

The state is the concatenation `[u₁; u₂; …; u_N]`. The parameters are a `CoupledParams`
holding the original parameters (`:p`), `:τ`, `:k`, `:adj`, etc.

## Example

```julia
using DynamicalSystems
function lorenz!(du, u, p, t)
    σ, ρ, β = p
    x, y, z = u
    du[1] = σ * (y - x)
    du[2] = x * (ρ - z) - y
    du[3] = x * y - β * z
    return nothing
end
ds = CoupledODEs(lorenz!, [1.0, 0.0, 0.0], [10.0, 28.0, 8/3])

cds = couple_systems( # 6-dimensional dynamical system
    ds, [1.0, 0.5], 1; k = 2.0,
    u0s = [[1.0, 0, 0], [-1.0, 2.0, 5.0]]
)

X, t = trajectory(cds, 50; Δt = 0.01)
x1, x2 = X[:, 1], X[:, 4]

set_parameter!(cds, :k, 0.5) # update coupling strength
set_parameter!(cds, :τ, [1.0, 0.8]) # update timescales
set_parameter!(cds, :p, [10.0, 35.0, 8/3]) # change parameters of the original system
````
"""
function couple_systems(
        ds::CoupledODEs, τ::AbstractVector{<:Real}, js;
        k = 1.0, adjacency = nothing, scale_coupling = false,
        u0s = [copy(current_state(ds)) for _ in 1:length(τ)], diffeq = ds.diffeq,
    )
    isinplace(ds) || error("Input system must be in-place (rule `f(du, u, p, t)`).")
    N = length(τ)
    N ≥ 2 || error("Need at least two systems to couple.")
    length(u0s) == N || error("Need one initial state per system.")
    D = dimension(ds)
    js_ = js isa Integer ? (Int(js),) : Tuple(Int.(js))
    all(j -> 1 ≤ j ≤ D, js_) || error("Coupled dimensions must lie in 1:$D.")

    adj = if isnothing(adjacency)
        [i == m ? 0.0 : 1.0 for i in 1:N, m in 1:N]
    else
        size(adjacency) == (N, N) || error("`adjacency` must be $N×$N.")
        Matrix{Float64}(adjacency)
    end

    T = eltype(current_state(ds))
    u0 = Vector{T}(reduce(vcat, Vector.(u0s)))

    p = CoupledParams(
        dynamic_rule(ds), current_parameters(ds), float.(collect(τ)),
        Float64(k), adj, js_, D, scale_coupling,
    )
    return CoupledODEs(coupled_rule!, u0, p; diffeq, t0 = initial_time(ds))
end