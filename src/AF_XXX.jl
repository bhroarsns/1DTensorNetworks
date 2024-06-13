include("./methods/itebd.jl")
include("./methods/idmrg.jl")

function hamiltonian(sitetype)
    modelname = "AF_XXX/$(replace(sitetype, '/'=>'_'))"
    i = addtags(siteind(sitetype), "i")
    j = addtags(siteind(sitetype), "j")
    hloc = op("Sx", i) * op("Sx", j) + op("Sy", i) * op("Sy", j) + op("Sz", i) * op("Sz", j)
    orginds = [i, j]
    return modelname, real(hloc), orginds
end

function executeTEBD(seed::Int, initΔτ::Float64, D::Int; initType="", mpslen=2, maxstep::Union{Int,Nothing}=nothing)
    modelname, hloc, orginds = hamiltonian("S=1/2")
    doiTEBD(modelname, hloc, orginds, "S=1/2", initΔτ, D, seed; mpslen, initType, maxstep)
    return nothing
end

function executeDMRG()
    modelname, hloc, orginds = hamiltonian("S=1/2")
    doiDMRG(modelname, hloc, orginds, "S=1/2", 16)
    return nothing
end