function tensorSV(mps::InfiniteMPS)
    svdict = Dict{String,Vector{Float64}}()
    for isite in 1:mps.length
        bl, _, bwl = bond(mps, isite - 1)
        bwr = bondWeight(mps, isite)
        st = siteTensor(mps, isite)
        si = siteInd(mps, isite)
        st2 = sqrt.(bwl) * st * sqrt.(bwr)
        for isi in eachindval(si)
            sarr = st2 * onehot(isi)
            _, _, _, spec = svd(sarr, [bl])
            svdict["$(sitename(mps, isite))$(isi[2])"] = spec.eigs
        end
    end
    return svdict
end

function compareSV(mps::InfiniteMPS, svdict::Dict{String,Vector{Float64}})
    diff = 0.0
    newdict = Dict{String,Vector{Float64}}()
    for isite in 1:mps.length
        bl, _, bwl = bond(mps, isite - 1)
        bwr = bondWeight(mps, isite)
        st = siteTensor(mps, isite)
        si = siteInd(mps, isite)
        st2 = sqrt.(bwl) * st * sqrt.(bwr)
        for isi in eachindval(si)
            smatr = st2 * onehot(isi)
            matrname = "$(sitename(mps, isite))$(isi[2])"
            _, _, _, spec = svd(smatr, [bl])
            newdict[matrname] = spec.eigs
            prelen = length(svdict[matrname])
            curlen = length(spec.eigs)
            truelen = max(prelen, curlen)
            curdiff = 0.0
            for isv in 1:truelen
                if isv <= prelen && isv <= curlen
                    curdiff += (svdict[matrname][isv] - spec.eigs[isv])^2
                else
                    if isv > prelen
                        curdiff += spec.eigs[isv]^2
                    else
                        curdiff += svdict[matrname][isv]^2
                    end
                end
            end
            diff += curdiff
        end
    end
    return diff, newdict
end