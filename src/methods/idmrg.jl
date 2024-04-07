include("../modules/DMRG.jl")
include("../modules/util.jl")

function logEigs(snapshotdir::String, curstep::Int, eigs::Vector{Float64}, irreps::Vector{Int})
    open("$(snapshotdir)/eigenvalues/$(curstep).dat", "w") do io
        for ieig in eachindex(eigs)
            println(io, "$(ieig), $(eigs[ieig]), $(irreps[ieig])")
        end
    end
end

function outputResult(resultdir::String, curstep::Int, gsEnergy::Float64, initialize=false)
    open("$(resultdir)/energy.dat", initialize ? "w" : "a") do io
        println(io, "$(curstep), $(gsEnergy/2/curstep)")
    end
end

function doiDMRG(
    # Hamiltonian
    modelname::String,
    hloc::ITensor,
    originalinds::Vector{Index{Int}},
    sitetype::String,
    # DMRG parameters
    D::Int;
    singlesite::Union{ITensor,Nothing}=nothing
)
    target = "$(modelname)/iDMRG/D=$(D)"
    resultdir, snapshotdir = setupDir(target)
    if !isdir("$(snapshotdir)/left")
        mkpath("$(snapshotdir)/left")
    end
    if !isdir("$(snapshotdir)/right")
        mkpath("$(snapshotdir)/right")
    end
    if !isdir("$(snapshotdir)/eigenvalues")
        mkpath("$(snapshotdir)/eigenvalues")
    end

    print("\r", 2)
    UL, UR, hL, hR, pleftbond, nleftsite, nleftbond, nrightbond, nrightsite, prightbond, eigs, irreps, gsEnergy = initDMRG(sitetype, hloc, originalinds, D; singlesite)
    outputResult(resultdir, 2, gsEnergy, true)
    snapshot("$(snapshotdir)/left/2.dat", UL, nleftsite, pleftbond, nleftbond)
    snapshot("$(snapshotdir)/right/2.dat", UR, nrightsite, prightbond, nrightbond)
    logEigs(snapshotdir, 2, eigs, irreps)

    for istep in 3:512
        print("\r", istep)
        pleftsite = nleftsite
        pleftbond = nleftbond
        prightbond = nrightbond
        prightsite = nrightsite
        UL, UR, hL, hR, nleftsite, nleftbond, nrightbond, nrightsite, eigs, irreps, gsEnergy = dmrgStep(istep, sitetype, hloc, originalinds, UL, UR, hL, hR, pleftbond, pleftsite, prightsite, prightbond, D; singlesite)
        outputResult(resultdir, istep, gsEnergy)
        snapshot("$(snapshotdir)/left/$(istep).dat", UL, nleftsite, pleftbond, nleftbond)
        snapshot("$(snapshotdir)/right/$(istep).dat", UR, nrightsite, prightbond, nrightbond)
        logEigs(snapshotdir, istep, eigs, irreps)
    end
    return nothing
end