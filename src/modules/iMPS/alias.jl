function getsitenum(mps::InfiniteMPS, sc::Char)
    sn = sc - 'A' + 1
    if 1 ≤ sn ≤ mps.length
        return sn
    else
        error("Invalid Site Character: \'$(sc)\'")
    end
end

function getbondnum(mps::InfiniteMPS, bn::String)
    leftsite = getsitenum(mps, bn[1])
    rightsite = getsitenum(mps, bn[2])
    if mod(leftsite + 1, 1:mps.length) == rightsite
        return leftsite
    else
        error("Invalid Bond Name: \"$(bondname)\"")
    end
end

function siteTensor(mps::InfiniteMPS, sc::Char)
    return siteTensor(mps, getsitenum(mps, sc))
end

function bondWeight(mps::InfiniteMPS, bn::String)
    return bondWeight(mps::InfiniteMPS, getbondnum(mps, bn))
end

function siteInd(mps::InfiniteMPS, sc::Char)
    return siteInd(mps, getsitenum(mps, sc))
end

function bondInds(mps::InfiniteMPS, sc::Char)
    return bondInds(mps, getsitenum(mps, sc))
end

function bond(mps::InfiniteMPS, bn::String)
    return bond(mps, getbondnum(mps, bn))
end

function contractKet(mps::InfiniteMPS, firstsite::Char, lastsite::Char; minketonly=false)
    return contractKet(mps, getsitenum(mps, firstsite), getsitenum(mps, lastsite); minketonly)
end

function transfermatrix(mps::InfiniteMPS, bn1::String, bn2::String=bn1)
    return transfermatrix(mps, getbondnum(mps, bn1), getbondnum(mps, bn2))
end

function environment(mps::InfiniteMPS, bn1::String, bn2::String=bn1)
    return environment(mps, getbondnum(mps, bn1), getbondnum(mps, bn2))
end

function densitymatrix(mps::InfiniteMPS, dmlen::Int, firstsite::Char; normalized=false)
    return densitymatrix(mps, dmlen, getsitenum(mps, firstsite); normalized)
end

function expectedvalue(mps::InfiniteMPS, op::ITensor, originalinds::Vector{Index{Int}}, firstsite::Char; normalized=false)
    return expectedvalue(mps, op, originalinds, getsitenum(mps, firstsite); normalized)
end

function canonicalize!(mps::InfiniteMPS, bn::String; fpcutoff::Float64=0.0, svcutoff::Float64=0.0)
    return canonicalize!(mps, getbondnum(mps, bn); fpcutoff, svcutoff)
end

function update!(mps::InfiniteMPS, gate::ITensor, originalinds::Vector{Index{Int}}, firstsite::Char; svcutoff::Float64=0.0, newbonddim::Union{Int,Nothing}=nothing)
    return update!(mps, gate, originalinds, getsitenum(mps, firstsite); svcutoff, newbonddim)
end