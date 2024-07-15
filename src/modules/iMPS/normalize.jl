function rescale!(mps::InfiniteMPS, λs::Vector{Float64})
    divider = sqrt(sum(abs.(λs)) / float(length(λs)))
    for ibond in eachindex(mps.bondWeights)
        bwnorm = norm(mps.bondWeights[ibond])
        divider /= bwnorm
        mps.bondWeights[ibond] /= bwnorm
    end
    divider ^= (1.0 / float(mps.length))
    for isite in eachindex(mps.siteTensors)
        mps.siteTensors[isite] /= divider
    end
    return nothing
end

function normalize!(mps::InfiniteMPS; opr::Dict{String,String}=Dict{String,String}(), plevel::UInt8=0b111)
    nopr = mergewith(*, opr, Dict("methodcall" => "normalize!,"))
    λs, curTrs = canonicalizeAll!(mps; opr=nopr, plevel)
    rescale!(mps, λs)
    return curTrs
end
