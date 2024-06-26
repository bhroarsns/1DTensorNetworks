function canonicalize!(mps::InfiniteMPS, bondnum::Int; opr::Dict{String,String}=Dict{String,String}(), svcutoff::Float64=0.0, parallel=true)
    nopr = mergewith(*, opr, Dict("methodcall" => "canonicalize!,"))
    mpslen = mps.length
    Dl, Ul, ll, Dr, Ur, lr, el, λ = environment(mps, bondnum; opr=nopr)

    bl, br, bw = bond(mps, bondnum)
    bn = bondname(mps, bondnum)

    Θ = sqrt.(Dl) * dag(Ul) * replaceinds(bw, [bl, br], [ll, lr]) * dag(Ur) * sqrt.(Dr)
    U, Σ, V = svd(Θ, el; cutoff=svcutoff, lefttags="Bond,$(bn),$(bn[1])", righttags="Bond,$(bn),$(bn[2])")

    let (errΘIO, prefix) = ssio(nopr, "errΘ")
        if !isnothing(errΘIO)
            println(errΘIO, prefix..., norm(Θ - U * Σ * V))
            flush(errΘIO)
        end
    end

    X = Ul * inv.(sqrt.(Dl)) * U
    Y = V * inv.(sqrt.(Dr)) * Ur

    if parallel
        return λ, Σ, replaceind(X, ll, bl), replaceind(Y, lr, br)
    end

    left = siteTensor(mps, bondnum) * replaceind(X, ll, bl)
    right = replaceind(Y, lr, br) * siteTensor(mps, bondnum + 1)

    mps.bondWeights[mod(bondnum, 1:mpslen)] = Σ
    mps.siteTensors[mod(bondnum, 1:mpslen)] = left
    mps.siteTensors[mod(bondnum + 1, 1:mpslen)] = right

    let (errCanon, prefix) = ssio(nopr, "errC")
        if !isnothing(errCanon)
            El, Er, llink, rlink = transfermatrix(mps, bondnum; opr=nopr)
            println(errCanon, prefix..., norm(El * δ(llink, llink') - λ * δ(uniqueinds(El, [llink, llink'])...)), ", ", norm(Er * δ(rlink, rlink') - λ * δ(uniqueinds(Er, [rlink, rlink'])...)))
            flush(errCanon)
        end
    end

    return λ
end

function canonicalizeAll!(mps::InfiniteMPS; opr::Dict{String,String}=Dict{String,String}(), svcutoff::Float64=0.0, parallel=true)
    nopr = mergewith(*, opr, Dict("methodcall" => "canonicalizeAll!,"))
    λs = Vector{Complex}(undef, mps.length)
    if parallel
        Σs = Vector{ITensor}(undef, mps.length)
        Xs = Vector{ITensor}(undef, mps.length)
        Ys = Vector{ITensor}(undef, mps.length)
        Threads.@threads for ibond in eachindex(λs)
            λ, Σ, X, Y = canonicalize!(mps, ibond; opr=merge(nopr, Dict("bond" => bondname(mps, ibond))), svcutoff, parallel)
            λs[ibond] = λ
            Σs[ibond] = Σ
            Xs[ibond] = X
            Ys[mod(ibond + 1, 1:mps.length)] = Y
        end
        for ibond in eachindex(Σs)
            mps.bondWeights[ibond] = Σs[ibond]
            mps.siteTensors[ibond] = Ys[ibond] * mps.siteTensors[ibond] * Xs[ibond]
        end
    else
        for ibond in eachindex(λs)
            λ = canonicalize!(mps, ibond; opr=merge(nopr, Dict("bond" => bondname(mps, ibond))), svcutoff, parallel)
            λs[ibond] = λ
        end
    end
    return λs
end

function normalize!(mps::InfiniteMPS; opr::Dict{String,String}=Dict{String,String}(), svcutoff::Float64=0.0, parallel=true)
    nopr = mergewith(*, opr, Dict("methodcall" => "normalize!,"))
    λs = canonicalizeAll!(mps; opr=nopr, svcutoff, parallel)
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