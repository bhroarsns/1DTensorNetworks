include("modules/iMPS.jl")

sitetype = "S=1/2"
bonddim1 = 16
bonddim2 = bonddim1
seed = 10

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


mps = randomMirrorInfiniteMPS(sitetype, bonddim1; seed)
