include("./methods/itebd.jl")
include("./methods/idmrg.jl")

function hamiltonian(;U::Float64, μ::Float64=-U/2.0)
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
    hloc += δ(j, j') * op("n↑ * n↓", i) * U / 2.0
    hloc += δ(i, i') * op("n↑ * n↓", j) * U / 2.0
    hloc += δ(j, j') * op("ntot", i) * -μ / 2.0
    hloc += δ(i, i') * op("ntot", j) * -μ / 2.0
    orginds = [i, j]
    return modelname, hloc, orginds
end

function executeTEBD(;U::Float64, μ::Float64=-U/2.0)
    modelname, hloc, orginds = hamiltonian(;U, μ)
    doTEBD(modelname, hloc, orginds, "Electron", [(0.1, 500)], 32, 10)
    return nothing
end

function executeDMRG(;U::Float64, μ::Float64=-U/2.0)
    modelname, hloc, orginds = hamiltonian(;U, μ)
    doiDMRG(modelname, hloc, orginds, "Electron", 16)
    return nothing
end