export tristable_fitzhugh_nagumo

function tristable_fitzhugh_nagumo(;
        β = 1, α = 1, e = 1.0,
        d = 5.0,
        a1 = 2.7778,
        a2 = 2.3333,
        a3 = 0.7619,
        a4 = 0.0847,
        c = 100.0,
        ε = 0.01,
        # have a look at Fig.2 of the paper for the regimes.
        # I use the tristable regime:
        a0 = -1,
        I = 0,
    )
    function tfn_rule!(du, u, p, t)
        v, w = u
        a0, ε, I, a1, a2, a3, a4, α, β, c, d, e = p
        du[1] = a0*v + a1*v^3 - a2*v^5 + a3*v^7 - a4*v^9 - α*w + I
        du[2] = ε*(β + c*v + d*v^2 - e*w)
        return nothing
    end
    # Note, I changed the order of parameters so that the most important three parameters are first
    p = [a0, ε, I, a1, a2, a3, a4, α, β, c, d, e]
    u0 = rand(2)
    diffeq = (abstol = 1e-6, reltol = 1e-6)
    return CoupledODEs(tfn_rule!, u0, p; diffeq)
end
