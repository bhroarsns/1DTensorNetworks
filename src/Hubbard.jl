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
    singlesite = nothing
    if U != 0.0
        singlesite = op("n↑ * n↓", i) * U
        singlesite += op("ntot", i) * -μ
    end
    orginds = [i, j]
    return modelname, hloc, orginds, singlesite
end

function executeTEBD(;U::Float64, μ::Float64=-U/2.0)
    modelname, hloc, orginds, singlesite = hamiltonian(;U, μ)
    doiTEBD(modelname, hloc, orginds, "Electron", [(0.1, 500)], 32, 10; singlesite)
    return nothing
end

function executeDMRG(;U::Float64, μ::Float64=-U/2.0)
    modelname, hloc, orginds, singlesite = hamiltonian(;U, μ)
    doiDMRG(modelname, hloc, orginds, "Electron", 16; singlesite)
    return nothing
end