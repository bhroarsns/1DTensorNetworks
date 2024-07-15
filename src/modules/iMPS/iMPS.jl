using Distributed
@everywhere using ITensors
@everywhere using LinearAlgebra
using Random
using HDF5

@everywhere mutable struct InfiniteMPS
    length::Int
    siteTensors::Vector{ITensor}
    bondWeights::Vector{ITensor}
end

@everywhere function siteTensor(mps::InfiniteMPS, sitenum::Int)
    return mps.siteTensors[mod(sitenum, 1:mps.length)]
end

@everywhere function bondWeight(mps::InfiniteMPS, bondnum::Int)
    return mps.bondWeights[mod(bondnum, 1:mps.length)]
end

@everywhere function sitename(mps::InfiniteMPS, sitenum::Int)
    return 'A' + mod(sitenum, 1:mps.length) - 1
end

@everywhere function bondname(mps::InfiniteMPS, bondnum::Int)
    return string(sitename(mps, bondnum), sitename(mps, bondnum + 1))
end

function correctTags!(mps::InfiniteMPS)
    for site in eachindex(mps.siteTensors)
        st = siteTensor(mps, site)
        bwl = bondWeight(mps, site - 1)
        bwr = bondWeight(mps, site)

        bl = commonind(st, bwl)
        si = uniqueind(st, bwl, bwr)
        br = commonind(st, bwr)
        oldinds = [bl, si, br]

        blnew = removetags(bl, "Site,TMLink")
        blnew = addtags(blnew, string("Bond,", sitename(mps, site), ", ", bondname(mps, site - 1)))
        sinew = removetags(si, "Bond,TMLink")
        sinew = addtags(sinew, string("Site,", sitename(mps, site)))
        brnew = removetags(br, "Site,TMLink")
        brnew = addtags(brnew, string("Bond,", sitename(mps, site), ", ", bondname(mps, site)))
        newinds = [blnew, sinew, brnew]

        replaceinds!(st, oldinds, newinds)
        replaceinds!(bwl, oldinds, newinds)
        replaceinds!(bwr, oldinds, newinds)
    end
end

function siteInd(mps::InfiniteMPS, sitenum::Int)
    return first(filter(hastags("Site"), inds(siteTensor(mps, sitenum))))
end

function bondInds(mps::InfiniteMPS, sitenum::Int)
    st = siteTensor(mps, sitenum)
    bwl = bondWeight(mps, sitenum - 1)
    bwr = bondWeight(mps, sitenum)
    return commonind(st, bwl), commonind(st, bwr)
end

@everywhere function bond(mps::InfiniteMPS, bondnum::Int)
    bw = bondWeight(mps, bondnum)
    stl = siteTensor(mps, bondnum)
    str = siteTensor(mps, bondnum + 1)
    return commonind(bw, stl), commonind(bw, str), bw
end

function InfiniteMPS(siteTensors::Vector{ITensor}, bondWeights::Vector{ITensor})
    if length(siteTensors) != length(bondWeights)
        error("Length of siteTensors and bondWeights must be the same.")
    end
    mps = InfiniteMPS(length(siteTensors), siteTensors, bondWeights)
    correctTags!(mps)
    return mps
end

function record(mps::InfiniteMPS, filename::String)
    f = h5open(filename, "w")
    for isite in eachindex(mps.siteTensors)
        write(f, "$(sitename(mps, isite))", mps.siteTensors[isite])
        write(f, "$(bondname(mps, isite))", mps.bondWeights[isite])
    end
    close(f)
    return nothing
end

include("initialState.jl")
@everywhere include("contractKet.jl")
@everywhere include("environment.jl")
@everywhere include("canonicalize.jl")
include("normalize.jl")
include("update.jl")
include("expectedValue.jl")
include("sv.jl")
include("correlation.jl")
include("alias.jl")

if !@isdefined ssio
    @everywhere function ssio(_::Dict{String,String}, _::String)
        return nothing, ""
    end
end
