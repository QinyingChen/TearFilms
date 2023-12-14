# uses Parameters pacakge
@with_kw struct PhysicalConstants
    p0 = 12.1e-6
    vw = 1.8e-5
    c0 = 300.0
    sigma_0 = 0.045
    mu = 1.3e-3
    d = 4.5e-6
    Df = 0.39e-9
    D0 = 1.6e-9
end

@with_kw struct PDEConstants
    ell::Float64
    eps::Float64
    Pc::Float64
    Pec::Float64
    invPec::Float64
end

function PDEConstants(p::PhysicalConstants, v_max::Real)
    @unpack p0, vw, c0, sigma_0, mu, d, D0 = p
    ell = (sigma_0 / mu / v_max)^(1 / 4) * d
    eps = d / ell
    Pc = (p0 * vw * c0) / v_max
    Pec = (v_max * ell) / (eps * D0)
    invPec = 1 / Pec
    return PDEConstants(ell, eps, Pc, Pec, invPec)
end
