function eachSiteCombination(mps::InfiniteMPS, t1::Int, t2::Int=t1)
    if t1 == t2
        si = siteInd(mps, t1)
        return Iterators.flatmap(is1 -> map(is2 -> (si => is1, si => is2), is1+1:dim(si)), eachval(si))
    else
        si1 = siteInd(mps, t1)
        si2 = siteInd(mps, t2)
        return Iterators.product(eachindval(si1), eachindval(si2))
    end
end

function divideBondWeights(mps::InfiniteMPS)
    siteTensors = deepcopy(mps.siteTensors)

    for it in 1:mps.length
        bl, br, bw = bond(mps, it)
        il = it
        ir = mod(it + 1, 1:mps.length)
        b = settags(bl, "Bond,$(bondname(mps,it))")
        siteTensors[il] = siteTensors[il] * replaceind(map(sqrt, bw), br, b)
        siteTensors[ir] = replaceind(map(sqrt, bw), bl, b) * siteTensors[ir]
    end

    return siteTensors
end

function commutator(sts::Vector{ITensor}, t1::Int, s1::Pair{<:Index,Int}, t2::Int, s2::Pair{<:Index,Int})
    mat1 = sts[t1] * onehot(s1)
    mat2 = sts[t2] * onehot(s2)
    return map(ci -> begin
        new1 = mat1 * replaceind(dag(mat1), ci, ci')
        new2 = mat2 * replaceind(dag(mat2), ci, ci')
        new12 = new1 * prime(new2)
        comm = new12 - new2 * prime(new1)
        replace(string(tags(removetags(ci, "Bond"))), "\"" => ""), norm(comm) / norm(new12)
    end, commoninds(mat1, mat2))
end

function correlation(mps::InfiniteMPS; opr::Dict{String,String}=Dict{String,String}())
    nopr = mergewith(*, opr, Dict("methodcall" => "correlation,"))
    sts = divideBondWeights(mps)
    for it in 1:mps.length, (is1, is2) in eachSiteCombination(mps, it)
        tname = sitename(mps, it)
        for (ibond, corr) in commutator(sts, it, is1, it, is2)
            let (corrIO, prefix) = ssio(merge(nopr, Dict("pair" => "$(tname)$(is1[2])-$(tname)$(is2[2])-$(ibond)")), "corr")
                if !isnothing(corrIO)
                    print(corrIO, prefix...)
                    println(corrIO, corr)
                    flush(corrIO)
                    close(corrIO)
                end
            end
        end
    end
    if mps.length == 2
        for (is1, is2) in eachSiteCombination(mps, 1, 2)
            for (ibond, corr) = commutator(sts, 1, is1, 2, is2)
                let (corrIO, prefix) = ssio(merge(nopr, Dict("pair" => "A$(is1[2])-B$(is2[2])-$(ibond)")), "corr")
                    if !isnothing(corrIO)
                        print(corrIO, prefix...)
                        println(corrIO, corr)
                        flush(corrIO)
                        close(corrIO)
                    end
                end
            end
        end
    else
        for it in 1:mps.length, (is1, is2) in eachSiteCombination(mps, it, it + 1)
            tname1 = sitename(mps, it)
            tname2 = sitename(mps, it + 1)
            for (is1, is2) in eachSiteCombination(mps, 1, 2)
                for (ibond, corr) = commutator(sts, it, is1, mod(it + 1, 1:mps.length), is2)
                    let (corrIO, prefix) = ssio(merge(nopr, Dict("pair" => "$(tname1)$(is1[2])-$(tname2)$(is2[2])-$(ibond)")), "corr")
                        if !isnothing(corrIO)
                            print(corrIO, prefix...)
                            println(corrIO, corr)
                            flush(corrIO)
                            close(corrIO)
                        end
                    end
                end
            end
        end
    end
    return nothing
end