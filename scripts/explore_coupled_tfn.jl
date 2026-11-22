using DrWatson
@quickactivate "MultistabilitySynchronization"

include(srcdir("couple_systems.jl"))
include(srcdir("predefined_systems.jl"))
include(srcdir("theme.jl"))

using GLMakie

ds = tristable_fitzhugh_nagumo()
u0s = [[0, y] for y in [0.1, 1.3, 2.0]]

# %% part one, interactive exploration of a single system
parameter_sliders = Dict(8 => 0.001:0.001:0.1)

fig, dsobs = interactive_trajectory_timeseries(ds, [1, 2], u0s; parameter_sliders, Δt = 0.05)

display(fig)

# %% part two, interactive exploration of a coupled system

cds = couple_systems( # 6-dimensional dynamical system
    ds, [1.02, 0.98], 1; k = 0.05, u0s = [u0s[1], u0s[2]],
)

u0sc = [[u; u .+ 1e-2randn(2)] for u in u0s]
parameter_sliders = Dict(:k => 0.001:0.001:0.5)

fig, dsobs = interactive_trajectory_timeseries(cds, [1, 3], u0sc;
    idxs = [1, 3], parameter_sliders, Δt = 0.05
)

# set_parameter!(cds, :k, 0.25) # change parameters of the original system

display(fig)
