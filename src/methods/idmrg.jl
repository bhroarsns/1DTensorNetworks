include("../modules/DMRG.jl")
include("../modules/util.jl")

function doiDMRG(
    # Hamiltonian
    modelname::String,
    hloc::ITensor,
    originalinds::Vector{Index{Int}},
    sitetype::String,
    # DMRG parameters
    D::Int,
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

    UL, UR, hL, hR, plb, nls, nlb, nrb, nrs, prb, eigs, irreps, gsEnergy = initDMRG(hloc, originalinds, sitetype, D)
    open("$(resultdir)/energy.dat", "w") do io
        println(io, "2, $(gsEnergy)")
    end
    snapshot("$(snapshotdir)/left/2.dat", UL, nls, plb, nlb)
    snapshot("$(snapshotdir)/right/2.dat", UR, nrs, prb, nrb)
    open("$(snapshotdir)/eigenvalues/2.dat", "w") do io
        for ieig in eachindex(eigs)
            println(io, "$(ieig), $(eigs[ieig]), $(irreps[ieig])")
        end
    end

    for istep in 3:512
        UL, UR, hL, hR, nlb, nls, nrs, nrb, eigs, irreps, gsEnergy = dmrgStep(
            istep,
            hloc,
            originalinds,
            sitetype,
            UL,
            UR,
            hL,
            hR,
            nlb,
            nls,
            nrs,
            nrb,
            D)

        open("$(resultdir)/energy.dat", "w") do io
            println(io, "$(istep), $(gsEnergy)")
        end
        snapshot("$(snapshotdir)/left/$(istep).dat", UL, nls, plb, nlb)
        snapshot("$(snapshotdir)/right/$(istep).dat", UR, nrs, prb, nrb)
        open("$(snapshotdir)/eigenvalues/$(istep).dat", "w") do io
            for ieig in eachindex(eigs)
                println(io, "$(ieig), $(eigs[ieig]), $(irreps[ieig])")
            end
        end
    end
end