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

    print("\r", 2)
    # plb - UL - nlb nrb - UR - prb
    #      nls             nrs
    UL, UR, hL, hR, plb, nls, nlb, nrb, nrs, prb, eigs, irreps, gsEnergy = initDMRG(hloc, originalinds, sitetype, D)
    open("$(resultdir)/energy.dat", "w") do io
        println(io, "2, $(gsEnergy/4)")
    end
    snapshot("$(snapshotdir)/left/2.dat", UL, nls, plb, nlb)
    snapshot("$(snapshotdir)/right/2.dat", UR, nrs, prb, nrb)
    open("$(snapshotdir)/eigenvalues/2.dat", "w") do io
        for ieig in eachindex(eigs)
            println(io, "$(ieig), $(eigs[ieig]), $(irreps[ieig])")
        end
    end

    for istep in 3:512
        print("\r", istep)
        #  UL - plb plb - UL - nlb nrb - UR - prb prb - UR
        # pls            nls             nrs            prs
        pls = nls
        plb = nlb
        prb = nrb
        prs = nrs
        UL, UR, hL, hR, nls, nlb, nrb, nrs, eigs, irreps, gsEnergy = dmrgStep(
            istep,
            hloc,
            originalinds,
            sitetype,
            UL,
            UR,
            hL,
            hR,
            plb,
            pls,
            prs,
            prb,
            D)

        open("$(resultdir)/energy.dat", "a") do io
            println(io, "$(istep), $(gsEnergy/2/istep)")
        end
        snapshot("$(snapshotdir)/left/$(istep).dat", UL, nls, plb, nlb)
        snapshot("$(snapshotdir)/right/$(istep).dat", UR, nrs, prb, nrb)
        open("$(snapshotdir)/eigenvalues/$(istep).dat", "w") do io
            for ieig in eachindex(eigs)
                println(io, "$(ieig), $(eigs[ieig]), $(irreps[ieig])")
            end
        end
    end
    return nothing
end