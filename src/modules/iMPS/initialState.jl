function randomInfiniteMPS(sitetype::String, bonddims::Vector{Int}; seed::Union{Int,Nothing}=nothing)
    mpslen = length(bonddims)
    siteTensors = Vector{ITensor}(undef, mpslen)
    bondWeights = Vector{ITensor}(undef, mpslen)
    bondindsl = Vector{Index{Int}}(undef, mpslen)
    bondindsr = Vector{Index{Int}}(undef, mpslen)
    for ibond in eachindex(bondWeights)
        bd = bonddims[ibond]
        lb = Index(bd)
        rb = Index(bd)
        bondWeights[ibond] = delta(lb, rb) / sqrt(bd)
        bondindsl[ibond] = lb
        bondindsr[ibond] = rb
    end
    if !isnothing(seed)
        Random.seed!(seed)
    end
    for isite in eachindex(siteTensors)
        si = siteind(sitetype)
        siteTensors[isite] = randomITensor(bondindsr[isite == 1 ? mpslen : isite - 1], si, bondindsl[isite])
    end
    return InfiniteMPS(siteTensors, bondWeights)
end

function randomInfiniteMPS(sitetype::String, bonddim::Int, mpslen::Int=2; seed::Union{Int,Nothing}=nothing)
    return randomInfiniteMPS(sitetype, ones(Int, mpslen) .* bonddim; seed)
end

function randomMirrorInfiniteMPS(sitetype::String, bonddim1::Int, bonddim2::Int=bonddim1; seed::Union{Int,Nothing}=nothing)
    sA = siteind(sitetype)
    sB = siteind(sitetype)
    bABA = Index(bonddim1)
    bABB = Index(bonddim1)
    bBAB = Index(bonddim2)
    bBAA = Index(bonddim1)

    if !isnothing(seed)
        Random.seed!(seed)
    end

    stA = randomITensor(bBAA, sA, bABA)
    stB = replaceinds(stA, [bBAA, sA, bABA], [bBAB, sB, bABB])
    bwAB = delta(bABA, bABB) / sqrt(bonddim1)
    bwBA = delta(bBAB, bBAA) / sqrt(bonddim2)

    return InfiniteMPS([stA, stB], [bwAB, bwBA])
end

function randomTIInfiniteMPS(sitetype::String, bonddim1::Int, bonddim2::Int=bonddim1; seed::Union{Int,Nothing}=nothing)
    sA = siteind(sitetype)
    sB = siteind(sitetype)
    bABA = Index(bonddim1)
    bABB = Index(bonddim1)
    bBAB = Index(bonddim2)
    bBAA = Index(bonddim1)

    if !isnothing(seed)
        Random.seed!(seed)
    end

    stA = randomITensor(bBAA, sA, bABA)
    stB = replaceinds(stA, [bBAA, sA, bABA], [bABB, sB, bBAB])
    bwAB = delta(bABA, bABB) / sqrt(bonddim1)
    bwBA = delta(bBAB, bBAA) / sqrt(bonddim2)

    return InfiniteMPS([stA, stB], [bwAB, bwBA])
end

function randomMirrorTIInfiniteMPS(sitetype::String, bonddim1::Int, bonddim2::Int=bonddim1; seed::Union{Int,Nothing}=nothing)
    sA = siteind(sitetype)
    sB = siteind(sitetype)
    bABA = Index(bonddim1)
    bABB = Index(bonddim1)
    bBAB = Index(bonddim2)
    bBAA = Index(bonddim1)

    if !isnothing(seed)
        Random.seed!(seed)
    end

    stA = randomITensor(bBAA, sA, bABA)
    stA = (stA + swapinds(stA, bBAA, bABA)) / 2.0
    stB = replaceinds(stA, [bBAA, sA, bABA], [bABB, sB, bBAB])
    bwAB = delta(bABA, bABB) / sqrt(bonddim1)
    bwBA = delta(bBAB, bBAA) / sqrt(bonddim2)

    return InfiniteMPS([stA, stB], [bwAB, bwBA])
end