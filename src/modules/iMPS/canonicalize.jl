function canonicalizer(mps::InfiniteMPS, bondnum::Int; opr::Dict{String,String}=Dict{String,String}(), plevel::UInt8=0b11)
    nopr = mergewith(*, opr, Dict("methodcall" => "canonicalizer,"))

    Dl, Ul, ll, Dr, Ur, lr, el, λ, curTr = environment(mps, bondnum; opr=nopr, plevel)
    bl, br, bw = bond(mps, bondnum)
    bn = bondname(mps, bondnum)

    Θ = sqrt.(Dl) * dag(Ul) * replaceinds(bw, [bl, br], [ll, lr]) * dag(Ur) * sqrt.(Dr)
    U, Σ, V = svd(Θ, el; lefttags="Bond,$(bn),$(bn[1])", righttags="Bond,$(bn),$(bn[2])")
    X = Ul * inv.(sqrt.(Dl)) * U
    Y = V * inv.(sqrt.(Dr)) * Ur

    let (errΘIO, prefix) = ssio(nopr, "errΘ")
        if !isnothing(errΘIO)
            println(errΘIO, prefix..., norm(Θ - U * Σ * V))
            flush(errΘIO)
        end
    end

    return λ, Σ, replaceind(X, ll, bl), replaceind(Y, lr, br), curTr
end

function canonicalize!(mps::InfiniteMPS, bondnum::Int; opr::Dict{String,String}=Dict{String,String}(), plevel::UInt8=0b11)
    nopr = mergewith(*, opr, Dict("methodcall" => "canonicalize!,"))
    mpslen = mps.length

    λ, Σ, X, Y, curTr = canonicalizer(mps, bondnum; opr=nopr, plevel)

    left = siteTensor(mps, bondnum) * X
    right = Y * siteTensor(mps, bondnum + 1)

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

    return λ, curTr
end

function canonicalizeAll!(mps::InfiniteMPS; opr::Dict{String,String}=Dict{String,String}(), plevel::UInt8=0b111)
    nopr = mergewith(*, opr, Dict("methodcall" => "canonicalizeAll!,"))
    if isodd(plevel)
        tasks = map(ibond -> remotecall(
                canonicalizer,
                1 + (ibond - 1) * (2^count_ones(plevel >> 1)),
                mps,
                ibond;
                opr=merge(nopr, Dict("bond" => bondname(mps, ibond))),
                plevel=plevel >> 1
            ), 1:2)
        results = fetch.(tasks)
        λs = getfield.(results, 1)
        Σs = getfield.(results, 2)
        Xs = getfield.(results, 3)
        Ys = getfield.(results, 4)
        curTrs = getfield.(results, 5)

        for ibond in eachindex(Σs)
            mps.bondWeights[ibond] = Σs[ibond]
            mps.siteTensors[ibond] = Ys[mod(ibond - 1, 1:mps.length)] * mps.siteTensors[ibond] * Xs[ibond]
        end
        return λs, curTrs
    else
        λs = zeros(mps.length)
        curTrs = zeros(mps.length)
        for ibond in eachindex(λs)
            λ, curTr = canonicalize!(mps, ibond; opr=merge(nopr, Dict("bond" => bondname(mps, ibond))), plevel=plevel >> 1)
            λs[ibond] = λ
            curTrs[ibond] = curTr
        end
        return λs, curTrs
    end
end
