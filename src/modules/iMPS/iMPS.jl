using ITensors
using LinearAlgebra
using Random

mutable struct InfiniteMPS
    length::Int
    siteTensors::Vector{ITensor}
    bondWeights::Vector{ITensor}
end

function siteTensor(mps::InfiniteMPS, sitenum::Int)
    return mps.siteTensors[mod(sitenum, 1:mps.length)]
end

function bondWeight(mps::InfiniteMPS, bondnum::Int)
    return mps.bondWeights[mod(bondnum, 1:mps.length)]
end

function sitename(mps::InfiniteMPS, sitenum::Int)
    return 'A' + mod(sitenum, 1:mps.length) - 1
end

function bondname(mps::InfiniteMPS, bondnum::Int)
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

function bond(mps::InfiniteMPS, bondnum::Int)
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

include("initialState.jl")
include("contractKet.jl")
include("environment.jl")
include("canonicalize.jl")
include("update.jl")
include("expectedValue.jl")
include("sv.jl")
include("correlation.jl")
include("alias.jl")

# function takeSnapshot(mps::InfiniteMPS; nopr::NamedTuple=(methodcall="",), ssio::Function=(_, _) -> (nothing, ""))
#     for isite in 1:mps.length
#         let (stIO, prefix) = ssio(nopr, "st:$(sitename(mps, isite))")
#             if !isnothing(stIO)
#                 st = mps.siteTensors[isite]
#                 si = siteInd(mps, isite)
#                 bl, br = bondInds(mps, isite)
#                 println(stIO, "# $(tags(si)), $(tags(bl)), $(tags(br)), real, imag, abs, angle")
#                 println(stIO, dim(si), ", ", dim(bl), ", ", dim(br))
#                 for isi in eachval(si)
#                     for ibl in eachval(bl)
#                         for ibr in eachval(br)
#                             entry = st[si=>isi, bl=>ibl, br=>ibr]
#                             println(stIO, isi, ", ", ibl, ", ", ibr, ", ", real(entry), ", ", imag(entry), ", ", abs(entry), ", ", angle(entry))
#                         end
#                     end
#                 end
#                 flush(stIO)
#             end
#         end
#     end

#     for ibond in 1:mps.length
#         let (bwIO, prefix) = ssio(nopr, "bw:$(bondname(mps, ibond))")
#             if !isnothing(bwIO)
#                 bl, _, bw = bond(mps, ibond)
#                 println(bwIO, "# $(bondname(mps, ibond)), real, imag, abs, angle")
#                 for ibl in eachval(bl)
#                     entry = bw[ibl, ibl]
#                     println(bwIO, ibl, ", ", real(entry), ", ", imag(entry), ", ", abs(entry), ", ", angle(entry))
#                 end
#                 flush(bwIO)
#             end
#         end
#     end
# end
