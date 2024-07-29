using Distributed
@everywhere function ssio(opr::Dict{String,String}, ssname::String)
    # ssname == "errtm" && return open("$(opr["snapshotdir"])/Err/$(opr["state"])/errtm.dat", "a"), (opr["step"], ", ")
    # ssname == "errv" && return open("$(opr["snapshotdir"])/Err/$(opr["state"])/errv_$(opr["side"]).dat", "a"), (opr["step"], ", ")
    # ssname == "errΘ" && return open("$(opr["snapshotdir"])/Err/$(opr["state"])/errΘ.dat", "a"), (opr["step"], ", ")
    # ssname == "errC" && return open("$(opr["snapshotdir"])/Err/$(opr["state"])/errC.dat", "a"), (opr["step"], ", ")
    ssname == "errU" && return open("$(opr["snapshotdir"])/Err/$(opr["state"])/errU.dat", "a"), (opr["step"], ", ", opr["fs"], ", ")
    ssname == "spec" && return open("$(opr["snapshotdir"])/Spec/$(opr["state"])/$(opr["bond"])_$(opr["side"])_$(opr["sector"]).dat", "a"), (opr["step"], ", ")
    # ssname == "degenFP" && return open("$(opr["snapshotdir"])/Step/$(opr["step"])/$(opr["state"])/$(opr["bond"])_$(opr["side"])_degenFP.dat", "w"), ""
    ssname == "corr" && return open("$(opr["snapshotdir"])/Corr/$(opr["pair"]).dat", "a"), (opr["step"], ", ")
    startswith(ssname, "uspec") && return open("$(opr["snapshotdir"])/Spec/$(opr["state"])/$(opr["fs"])_$(opr["bond"]).dat", "a"), (opr["step"], ", ")
    # if (opr["methodcall"] == "normalize!,") || (opr["methodcall"] == "update!,")
    #     startswith(ssname, "st") && return open("$(opr["snapshotdir"])/Step/$(opr["step"])/$(opr["state"])/$(ssname[4:end])", "w"), ""
    #     startswith(ssname, "bw") && return open("$(opr["snapshotdir"])/Step/$(opr["step"])/$(opr["state"])/$(ssname[4:end])", "w"), ""
    # end
    return nothing, ""
end

include("../modules/iMPS/iMPS.jl")
include("../modules/util.jl")
using Printf

function measurement(resultdir::String, mps::InfiniteMPS, hloc::ITensor, originalinds::Vector{Index{Int}}, istep::Int, β::Float64; singlesite::Union{ITensor,Nothing}=nothing, obs::Union{Vector{Tuple{ITensor,Vector{Index{Int}}}},Nothing}=nothing)
    bondEs = real.(expectedvalues(mps, hloc, originalinds; normalized=true))
    aveE = sum(bondEs) / length(bondEs)
    output = @sprintf("%03d, %+.16e, %+.16e, %s", istep, β, aveE, join(map(ev -> @sprintf("%+.16e", ev), bondEs), ", "))
    if !isnothing(singlesite)
        ssEs = real.(expectedvalues(mps, singlesite, [originalinds[begin]]; normalized=true))
        aveE += sum(ssEs) / length(ssEs)
        output = @sprintf("%03d, %+.16e, %+.16e, %s, %s", istep, β, aveE, join(map(ev -> @sprintf("%+.16e", ev), bondEs), ", "), join(map(ev -> @sprintf("%+.16e", ev), ssEs), ", "))
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
    return aveE
end

function printSV(snapshotdir::String, svdict::Dict{String,Vector{Float64}})
    for k in keys(svdict)
        open("$(snapshotdir)/$(k).dat", "a") do io
            println(io, join(svdict[k], ", "))
        end
    end
    return nothing
end

function diffSV(snapshotdir::String, mps::InfiniteMPS, svdict::Dict{String,Vector{Float64}})
    diff, newdict = compareSV(mps, svdict)
    printSV(snapshotdir, newdict)
    return diff, newdict
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
    # optional Hamiltonian parameters
    singlesite::Union{ITensor,Nothing}=nothing,
    obs::Union{Vector{Tuple{ITensor,Vector{Index{Int}}}},Nothing}=nothing,
    # optional TEBD parameters
    fullspec=true,
    maxstep::Union{Int,Nothing}=nothing,
    haltthres::Float64=1.0e-10,
    numΔτ::Int=20,
    eneThreas::Float64=1.0e-6,
    # optional MPS parameters
    mpslen::Int=length(originalinds),
    initType="",
    # optional computational parameters
    verbose=true,
    plevel::UInt8=0b111,
    recordInterval::Union{Int,Nothing}=nothing,
    date::String=string(Date(Dates.now())),
)
    # directory setup
    target = "$(modelname)/iTEBD/mpslen=$(mpslen)/D=$(D)/seed=$(seed)/initΔτ=$(initΔτ)" * (!isempty(initType) ? "/$(replace(initType, '/' => '-'))" : "")
    resultdir, snapshotdir = setupDir(date, target)
    mkpathINE("$(snapshotdir)/Corr")
    mkpathINE("$(snapshotdir)/Spec/BSU")
    mkpathINE("$(snapshotdir)/Err/BSU")
    mkpathINE("$(snapshotdir)/Spec/FUN")
    mkpathINE("$(snapshotdir)/Err/FUN")
    if !isnothing(singlesite)
        mkpathINE("$(snapshotdir)/Spec/BSN")
        mkpathINE("$(snapshotdir)/Err/BSN")
        mkpathINE("$(snapshotdir)/Spec/FUU")
        mkpathINE("$(snapshotdir)/Err/FUU")
    end
    mkpathINE("$(snapshotdir)/Step/0/FUN")
    touch("$(resultdir)/convergence.dat")

    logtime = open("./logTime.txt", "w")
    println(logtime, target); flush(logtime)

    # parameter initialization
    β = 0.0
    Δτ = initΔτ
    totsteps = 0
    opr = Dict("snapshotdir" => snapshotdir, "step" => string(0), "methodcall" => "")

    # mps initialization
    mps = if initType == "Mirror"
        randomMirrorInfiniteMPS(sitetype, D; seed)
    elseif initType == "TI"
        randomTIInfiniteMPS(sitetype, D; seed)
    elseif initType == "MirrorTI"
        randomMirrorTIInfiniteMPS(sitetype, D; seed)
    elseif initType == "SII"
        randomSIIInfiniteMPS(sitetype, D, mpslen; seed)
    elseif initType == "PHS"
        randomPHSInfiniteMPS(sitetype, D, mpslen; seed)
    elseif initType == "PHESII"
        randomPHESIIInfiniteMPS(sitetype, D, mpslen; seed)
    else
        if isempty(initType)
            randomInfiniteMPS(sitetype, D, mpslen; seed)
        else
            fromHDF5(initType)
        end
    end

    # init state measurements
    curTrs = normalize!(mps; opr=merge(opr, Dict("state" => "FUN")), plevel)
    prevsv = tensorSV(mps)
    printSV(snapshotdir, prevsv)
    prevE = measurement(resultdir, mps, hloc, originalinds, totsteps, β; singlesite, obs)

    for _ in 1:numΔτ
        gate = exp(-Δτ * hloc)
        sgate = isnothing(singlesite) ? nothing : exp(-Δτ * singlesite)
        istep = 0
        diff = Inf
        curE = prevE
        diffE = Inf

        while diff > haltthres || diffE > eneThreas
            # step number incrementation
            istep += 1
            curstep = totsteps + istep
            if !isnothing(maxstep) && curstep > maxstep
                open("$(resultdir)/convergence.dat", "a") do io
                    println(io, Δτ, ", ", curE, ", ", totsteps)
                end
                verbose && println("\n", "Interrupted due to maxstep($(maxstep))")
                @goto interuption
            end
            verbose && print("\r", curstep, ", ", Δτ)
            opr["step"] = string(curstep)

            # bond hamiltonian update
            mkpathINE("$(snapshotdir)/Step/$(curstep)/BSU")
            update!(mps, gate, originalinds; opr=merge(opr, Dict("state" => "BSU")))

            # single site gate update (if exists)
            if !isnothing(sgate)
                mkpathINE("$(snapshotdir)/Step/$(curstep)/BSN")
                print(logtime, (@elapsed normalize!(mps; opr=merge(opr, Dict("state" => "BSN")), plevel)), ", "); flush(logtime)
                mkpathINE("$(snapshotdir)/Step/$(curstep)/FUU")
                update!(mps, sgate, [originalinds[begin]]; opr=merge(opr, Dict("state" => "FUU")))
            end

            # canonicalization & normalization
            mkpathINE("$(snapshotdir)/Step/$(curstep)/FUN")
            print(logtime, (@elapsed normalize!(mps; opr=merge(opr, Dict("state" => "FUN")), plevel)), ", "); flush(logtime)

            # measurements
            correlation(mps; opr=merge(opr, Dict("state" => "FUN")))
            diff, prevsv = diffSV(snapshotdir, mps, prevsv)
            tmpE = measurement(resultdir, mps, hloc, originalinds, curstep, β + Δτ * istep; singlesite, obs)

            # monitoring
            println(logtime, @sprintf("%.16e, %.16e", diff, diffE)); flush(logtime)
            verbose && @printf ", total: %.16e, %.16e" diff diffE
            if !isnothing(recordInterval) && (curstep % recordInterval == 0)
                record(mps, "./$(snapshotdir)/Step/$(curstep)/mps.h5")
            end

            diffE = abs(tmpE - curE)
            curE = tmpE
        end
        open("$(resultdir)/convergence.dat", "a") do io
            println(io, Δτ, ", ", curE, ", ", totsteps+istep)
        end
        if abs(curE - prevE) / abs(prevE) < eneThreas
            verbose && println("\n", "Interrupted due to energy convergence(prevE=$(prevE), curE=$(curE))")
            @goto interuption
        end

        β += Δτ * istep
        totsteps += istep
        Δτ /= 2.0
        prevE = curE
        verbose && println("")
    end
    verbose && println("\n", "Interrupted due to numΔτ($numΔτ)")
    @label interuption

    println("\r", modelname, ": finished")
    close(logtime)
    record(mps, "./$(snapshotdir)/mps.h5")

    return nothing
end