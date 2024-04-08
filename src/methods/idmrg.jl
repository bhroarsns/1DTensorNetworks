include("../modules/DMRG.jl")
include("../modules/util.jl")

function logEigs(filename::String, eigs::Vector{Float64}, degen::Vector{Int})
    open(filename, "w") do io
        for ieig in eachindex(eigs)
            println(io, "$(ieig), $(eigs[ieig]), $(degen[ieig])")
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
    mkpathINE("$(snapshotdir)/left")
    mkpathINE("$(snapshotdir)/right")
    mkpathINE("$(snapshotdir)/eigenvalues")

    print("\r", 2)
    UL, UR, hL, hR, plb, nls, nlb, nrb, nrs, prb, eigs, degen, gsEnergy = initDMRG(sitetype, hloc, originalinds, D; singlesite)
    outputResult(resultdir, 2, gsEnergy, true)
    snapshot("$(snapshotdir)/left/2.dat", UL, nls, plb, nlb)
    snapshot("$(snapshotdir)/right/2.dat", UR, nrs, prb, nrb)
    logEigs("$(snapshotdir)/eigenvalues/2.dat", eigs, degen)

    for istep in 3:512
        print("\r", istep)
        pls = nls
        plb = nlb
        prb = nrb
        prs = nrs
        UL, UR, hL, hR, nls, nlb, nrb, nrs, eigs, degen, gsEnergy = dmrgStep(istep, sitetype, hloc, originalinds, UL, UR, hL, hR, plb, pls, prs, prb, D; singlesite)
        outputResult(resultdir, istep, gsEnergy)
        snapshot("$(snapshotdir)/left/$(istep).dat", UL, nls, plb, nlb)
        snapshot("$(snapshotdir)/right/$(istep).dat", UR, nrs, prb, nrb)
        logEigs("$(snapshotdir)/eigenvalues/$(istep).dat", eigs, degen)
    end
    return nothing
end