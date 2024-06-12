include("./methods/itebd.jl")
include("./methods/idmrg.jl")

function hamiltonian(; U::Float64, μ::Float64=U / 2.0)
    modelname = "Hubbard/U=$(U)"
    i = addtags(siteind("Electron"), "i")
    j = addtags(siteind("Electron"), "j")

    # site order: left to right 1, 2, 3, ...
    # op. order: younger first ... c†₂ c†₁|0>
    # spin order: up - down c†↑ c†↓ |0>
    # ==> j↑ j↓ i↑ i↓
    hloc = op("c†↑ * F↓", i) * op("c↑", j) * -1.0
    hloc += op("c†↓ * F↑", i) * op("c↓", j) * -1.0
    hloc += op("c↑ * F↓", i) * op("c†↑", j) * -1.0
    hloc += op("c↓ * F↑", i) * op("c†↓", j) * -1.0
    hbond = hloc
    singlesite = nothing
    if U != 0.0
        singlesite = op("n↑ * n↓", i) * U
        singlesite += op("ntot", i) * -μ
        singlesite += δ(i, i') * (U / 2.0)
        hbond += singlesite * δ(j, j') * 0.5
        hbond += δ(i, i') * replaceinds(singlesite, [i, i'], [j, j']) * 0.5
    end
    orginds = [i, j]
    return modelname, hloc, orginds, singlesite, hbond
end

function executeTEBD(seed::Int, initΔτ::Float64, D::Int; U::Float64, μ::Float64=U / 2.0, symmetry="")
    modelname, hloc, orginds, singlesite, hbond = hamiltonian(; U, μ)
    doiTEBD(
        modelname,
        hloc,
        orginds,
        "Electron",
        initΔτ,
        D,
        seed;
        singlesite,
        obs=[
            (op("ntot", orginds[begin]), [orginds[begin]]),
            (op("n↑", orginds[begin]), [orginds[begin]]),
            (op("n↓", orginds[begin]), [orginds[begin]]),
        ],
        symmetry)
    return nothing
end

function executeDMRG(; U::Float64, μ::Float64=U / 2.0)
    modelname, hloc, orginds, singlesite = hamiltonian(; U, μ)
    doiDMRG(modelname, hloc, orginds, "Electron", 8; singlesite)
    return nothing
end