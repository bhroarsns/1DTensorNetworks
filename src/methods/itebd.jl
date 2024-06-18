using Pkg; Pkg.activate(".")

include("../modules/iMPS.jl")
include("../modules/util.jl")
using Printf
using HDF5

function measurement(resultdir::String, mps::InfiniteMPS, hloc::ITensor, originalinds::Vector{Index{Int}}, istep::Int, β::Float64; singlesite::Union{ITensor,Nothing}=nothing, obs::Union{Vector{Tuple{ITensor,Vector{Index{Int}}}},Nothing}=nothing)
    evs = real.(expectedvalues(mps, hloc, originalinds; normalized=true))
    aveev = sum(evs) / length(evs)
    output = @sprintf("%03d, %+.16e, %+.16e, %s", istep, β, aveev, join(map(ev -> @sprintf("%+.16e", ev), evs), ", "))
    if !isnothing(singlesite)
        ssevs = real.(expectedvalues(mps, singlesite, [originalinds[begin]]; normalized=true))
        aveev += sum(ssevs) / length(ssevs)
        output = @sprintf("%03d, %+.16e, %+.16e, %s, %s", istep, β, aveev, join(map(ev -> @sprintf("%+.16e", ev), evs), ", "), join(map(ev -> @sprintf("%+.16e", ev), ssevs), ", "))
    end
    if !isnothing(obs)
        for (ob, orginds) in obs
            obevs = real.(expectedvalues(mps, ob, orginds; normalized=true))
            output *= @sprintf(", %+.16e, %s", sum(obevs) / length(obevs), join(map(ev -> @sprintf("%+.16e", ev), obevs), ", "))
        end
    end
    open("$(resultdir)/energy.dat", "a") do io
        println(io, output)
    end
end

function genssio(snapshotdir::String)
    return (opr::NamedTuple, ssname::String) -> begin
        mkpathINE("$(snapshotdir)/Step/$(opr.step)/$(opr.state)")
        mkpathINE("$(snapshotdir)/Spec/$(opr.state)")
        mkpathINE("$(snapshotdir)/Err/$(opr.state)")
        mkpathINE("$(snapshotdir)/Corr")
        # ssname == "errtm" && return open("$(snapshotdir)/Err/$(opr.state)/errtm.dat", "a"), (opr.step, ", ")
        ssname == "sspec" && return open("$(snapshotdir)/Spec/$(opr.state)/$(opr.bond)_$(opr.side)_sym.dat", "a"), (opr.step, ", ")
        ssname == "aspec" && return open("$(snapshotdir)/Spec/$(opr.state)/$(opr.bond)_$(opr.side)_asym.dat", "a"), (opr.step, ", ")
        ssname == "totspec" && return open("$(snapshotdir)/Spec/$(opr.state)/$(opr.bond)_$(opr.side)_tot.dat", "a"), (opr.step, ", ")
        # ssname == "errσ" && return open("$(snapshotdir)/Err/$(opr.state)/errσ.dat", "a"), (opr.step, ", ")
        # ssname == "errμ" && return open("$(snapshotdir)/Err/$(opr.state)/errμ.dat", "a"), (opr.step, ", ")
        ssname == "degenFPl" && return open("$(snapshotdir)/Step/$(opr.step)/$(opr.state)/$(opr.bond)_left_degenFP.dat", "w"), ""
        ssname == "degenFPr" && return open("$(snapshotdir)/Step/$(opr.step)/$(opr.state)/$(opr.bond)_right_degenFP.dat", "w"), ""
        # ssname == "errΘ" && return open("$(snapshotdir)/Err/$(opr.state)/errΘ.dat", "a"), (opr.step, ", ")
        # ssname == "errC" && return open("$(snapshotdir)/Err/$(opr.state)/errC.dat", "a"), (opr.step, ", ")
        ssname == "errU" && return open("$(snapshotdir)/Err/$(opr.state)/errU.dat", "a"), (opr.step, ", ", opr.fs, ", ")
        ssname == "corr" && return open("$(snapshotdir)/Corr/$(opr.pair).dat", "a"), (opr.step, ", ")
        startswith(ssname, "uspec") && return open("$(snapshotdir)/Spec/$(opr.state)/$(opr.fs)_$(opr.bond).dat", "a"), (opr.step, ", ")
        # if (opr.methodcall == "normalize!,") || (opr.methodcall == "update!,")
        #     startswith(ssname, "st") && return open("$(snapshotdir)/Step/$(opr.step)/$(opr.state)/$(ssname[4:end])", "w"), ""
        #     startswith(ssname, "bw") && return open("$(snapshotdir)/Step/$(opr.step)/$(opr.state)/$(ssname[4:end])", "w"), ""
        # end
        return nothing, nothing
    end
end

function printSV(snapshotdir::String, svdict::Dict{String, Vector{Float64}})
    for k in keys(svdict)
        open("$(snapshotdir)/$(k).dat", "a") do io
            println(io, join(svdict[k], ", "))
        end
    end
end

function doiTEBD(
    # Hamiltonian
    modelname::String,
    hloc::ITensor,
    originalinds::Vector{Index{Int}},
    sitetype::String,
    # TEBD parameters
    initΔτ::Float64,
    # MPS parameters
    D::Int,
    seed::Int;
    mpslen::Int=length(originalinds),
    singlesite::Union{ITensor,Nothing}=nothing,
    obs::Union{Vector{Tuple{ITensor,Vector{Index{Int}}}},Nothing}=nothing,
    initType="",
    maxstep::Union{Int,Nothing}=nothing
)
    target = "$(modelname)/iTEBD/mpslen=$(mpslen)/D=$(D)/seed=$(seed)/initΔτ=$(initΔτ)"*(!isempty(initType) ? "/$(initType)" : "")
    resultdir, snapshotdir = setupDir(target)
    open("$(resultdir)/energy.dat", "w") do io
        println(io, "# D=$(D), seed=$(seed)")
    end

    β = 0.0
    totsteps = 0
    mps = if initType == "Mirror"
        randomMirrorInfiniteMPS(sitetype, D; seed)
    elseif initType == "TI"
        randomTIInfiniteMPS(sitetype, D; seed)
    else
        randomInfiniteMPS(sitetype, D, mpslen; seed)
    end

    normalize!(mps; opr=(step=0, methodcall="", state="FUN"), ssio=genssio(snapshotdir))
    prevsv = tensorSV(mps)
    printSV(snapshotdir, prevsv)
    measurement(resultdir, mps, hloc, originalinds, totsteps, β; singlesite, obs)

    Δτ = initΔτ
    while initΔτ / Δτ < 1000.0
        gate = exp(-Δτ * hloc)
        sgate = isnothing(singlesite) ? nothing : exp(-Δτ * singlesite)
        istep = 0
        diff = Inf
        while diff > 1.0e-10
            istep += 1
            curstep = totsteps + istep
            if !isnothing(maxstep) && curstep > maxstep
                println("\r", modelname, ": interrupted")
                @goto end_of_loop
            end
            print("\r", curstep, ", ", Δτ)
            update!(mps, gate, originalinds; opr=(step=curstep, methodcall="", state="BSU"), ssio=genssio(snapshotdir))
            if !isnothing(sgate)
                normalize!(mps; opr=(step=curstep, methodcall="", state="BSN"), ssio=genssio(snapshotdir))
                update!(mps, sgate, [originalinds[begin]]; opr=(step=curstep, methodcall="", state="FUU"), ssio=genssio(snapshotdir))
            end
            normalize!(mps; opr=(step=curstep, methodcall="", state="FUN"), ssio=genssio(snapshotdir))
            correlation(mps; opr=(step=curstep, methodcall="", state="FUN"), ssio=genssio(snapshotdir))
            diff, prevsv = compareSV(mps, prevsv)
            printSV(snapshotdir, prevsv)
            @printf ", total: %.16e" diff
            measurement(resultdir, mps, hloc, originalinds, curstep, β + Δτ * istep; singlesite, obs)
        end
        β += Δτ * istep
        totsteps += istep
        Δτ /= 2.0
        println("")
    end
    @label end_of_loop
    println("\r", modelname, ": finished")

    f = h5open("./$(snapshotdir)/mps.h5", "w")
    for isite in eachindex(mps.siteTensors)
        write(f, "$(sitename(mps, isite))", mps.siteTensors[isite])
        write(f, "$(bondname(mps, isite))", mps.bondWeights[isite])
    end
    close(f)

    El, _, llink, _ = transfermatrix(mps, 1)
    leftInds = uniqueinds(El, [llink, llink'])
    P, Pinv, Q, Qinv, indsym, indasym = symprojector(llink, llink')
    symtm = P * El * replaceinds(Pinv, [llink, llink'], leftInds)
    asymtm = Q * El * replaceinds(Qinv, [llink, llink'], leftInds)
    open("$(snapshotdir)/tm.dat", "w") do io
        println(io, dim(llink), ", ", dim(llink'), ", ", dim(leftInds[1]), ", ", dim(leftInds[2]))
        for il1 in eachval(llink)
            for il2 in eachval(llink')
                for ir1 in eachval(leftInds[1])
                    for ir2 in eachval(leftInds[2])
                        entry = El[llink=>il1, llink'=>il2, leftInds[1]=>ir1, leftInds[2]=>ir2]
                        println(io, il1, ", ", il2, ", ", ir1, ", ", ir2, ", ", il1 + dim(llink) * (il2 - 1), ", ", ir1 + dim(leftInds[1]) * (ir2 - 1), ", ", real(entry), ", ", imag(entry), ", ", abs(entry), ", ", angle(entry))
                    end
                end
            end
        end
    end
    open("$(snapshotdir)/symtm.dat", "w") do io
        println(io, dim(indsym), ", ", dim(indsym'))
        for is1 in eachval(indsym)
            for is2 in eachval(indsym')
                entry = symtm[indsym=>is1, indsym'=>is2]
                println(io, is1, ", ", is2, ", ", real(entry), ", ", imag(entry), ", ", abs(entry), ", ", angle(entry))
            end
        end
    end
    open("$(snapshotdir)/asymtm.dat", "w") do io
        println(io, dim(indasym), ", ", dim(indasym'))
        for ia1 in eachval(indasym)
            for ia2 in eachval(indasym')
                entry = asymtm[indasym=>ia1, indasym'=>ia2]
                println(io, ia1, ", ", ia2, ", ", real(entry), ", ", imag(entry), ", ", abs(entry), ", ", angle(entry))
            end
        end
    end

    return nothing
end