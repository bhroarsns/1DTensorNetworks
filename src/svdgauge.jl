using ITensors
using HDF5

function snapshot(targetDir::String, name::String, tensor::ITensor, s::Index, bondAB::Index, bondBA::Index)
    for is in eachval(s)
        open("./snapshots/$(targetDir)/$(name)_$(is).dat", "w") do io
            println(io, dim(bondAB), " ", dim(bondBA))
            for ib in eachindval(bondAB, bondBA)
                entry = tensor[s=>is, ib...]
                println(io, join(getfield.(ib, 2), ", "), ", ", abs(entry), ", ", angle(entry), ", ", real(entry), ", ", imag(entry))
            end
        end
    end
end

function logTM(targetDir::String, name::String, tensorA::ITensor, tensorB::ITensor, bondBA::Index, sectoredInd::Vector{NTuple{4, Int}})
    linkBA = replacetags(bondBA, "Perm"=>"TMLink")
    ket = replaceind(tensorA, bondBA, linkBA) * tensorB
    tm = ket * noprime(prime(dag(ket)); tags="Site")
    combinedSectorInds = map(z -> z[2], sort!(collect(Iterators.flatmap(x -> map(y -> ((abs(x[2][1] - y[2][1]), abs(x[2][2] - y[2][2]), x[2][1], y[2][1], x[2][2], y[2][2], x[2][3], y[2][3], x[2][4], y[2][4]), (x[1], y[1])), Iterators.enumerate(sectoredInd)), Iterators.enumerate(sectoredInd)))))
    open("./snapshots/$(targetDir)/$(name)_tm.dat", "w") do io
        for lcsi in Iterators.enumerate(combinedSectorInds)
            for rcsi in Iterators.enumerate(combinedSectorInds)
                entry = tm[linkBA => lcsi[2][1], linkBA' => lcsi[2][2], bondBA => rcsi[2][1], bondBA' => rcsi[2][2]]
                println(io, lcsi[1], ", ", rcsi[1], ", ", abs(entry), ", ", angle(entry), ", ", real(entry), ", ", imag(entry))
            end
        end
    end
end

function divideBondWeights(stA::ITensor, stB::ITensor, bwAB::ITensor, bwBA::ITensor)
    A_BA, siteA, A_AB = inds(stA)
    B_AB, siteB, B_BA = inds(stB)

    bondAB = removetags(sim(A_AB), "A")
    bondBA = removetags(sim(A_BA), "A")

    A = replaceinds(sqrt.(bwBA) * stA * sqrt.(bwAB), [B_AB, B_BA], [bondAB, bondBA])
    B = replaceinds(sqrt.(bwAB) * stB * sqrt.(bwBA), [A_AB, A_BA], [bondAB, bondBA])

    return A, B, siteA, siteB, bondAB, bondBA
end

function svdGauge(A::ITensor, B::ITensor, isite::Int, siteA::Index{Int}, bondBA::Index{Int})
    matA = A * onehot(siteA => isite)
    U, S, V, _, linkBA, linkAB = svd(matA, [bondBA]; lefttags="Link,BA", righttags="Link,AB")
    A_sg = dag(U) * A * dag(V)
    B_sg = V * B * U

    return A_sg, B_sg, linkBA, linkAB
end

function svdDiagGauge(A_sg::ITensor, B_sg::ITensor, linkBA::Index{Int}, linkAB::Index{Int}, isite::Int, siteA::Index{Int}, siteB::Index{Int})
    matB = B_sg * onehot(siteB => dim(siteA) - isite + 1)
    D, P, _, eigenAB, eigenBA = eigen(matB, [linkAB], [linkBA]; lefttags="Eigen,AB", righttags="Eigen,BA")
    Pinv = ITensor(inv(matrix(P)), eigenAB, linkAB)

    A_sdg = replaceinds(Pinv, [eigenAB, linkAB], [eigenBA, linkBA]) * A_sg * replaceinds(P, [eigenBA, linkBA], [eigenAB, linkAB])
    B_sdg = Pinv * B_sg * P
    return A_sdg, B_sdg, eigenBA, eigenAB
end

    
function getConnection(matB::ITensor, D::Int)
    strM = abs.(storage(matB))
    logstrM = sort!(log10.(strM); rev=true)
    difflsM = map(i -> logstrM[i] - logstrM[i+1], 1:length(strM)-1)
    numNZ = findfirst(isequal(floor(maximum(difflsM))), floor.(difflsM))
    sp = map(iall -> ((iall - 1) % D + 1, (iall - 1) ÷ D + 1), sortperm(strM; rev=true)[1:numNZ])
    return sp
end

function getSectorPerm(B_sdg::ITensor, isite::Int, siteB::Index{Int}, eigenBA::Index{Int})
    D = dim(eigenBA)
    jsite = mod(isite+2, 1:4)
    matB1 = B_sdg * onehot(siteB => isite)
    matB2 = B_sdg * onehot(siteB => jsite)
    sp1 = getConnection(matB1, D)
    sp2 = getConnection(matB2, D)
    
    sectors = zeros(Int, D)
    numsecs = 0
    for iall in sp1
        i1 = iall[1]
        i2 = iall[2]
        if sectors[i1] == 0
            if sectors[i2] == 0
                numsecs += 1
                sectors[i1] = numsecs
                sectors[i2] = numsecs
            else
                sectors[i1] = sectors[i2]
            end
        else
            if sectors[i2] == 0
                sectors[i2] = sectors[i1]
            elseif sectors[i1] != sectors[i2]
                if sectors[i1] < sectors[i2]
                    sectors = map(x -> x > sectors[i2] ? x - 1 : x == sectors[i2] ? sectors[i1] : x, sectors)
                    numsecs -= 1
                else
                    sectors = map(x -> x > sectors[i1] ? x - 1 : x == sectors[i1] ? sectors[i2] : x, sectors)
                    numsecs -= 1
                end
            end
        end
    end
    # sectorNames = sort(unique(sectors))
    # numsecs = length(sectorNames)
    # replace!(sectors, map(x -> x[2] == 0 ? 0 => numsecs : x[2] => x[1], Iterators.enumerate(sectorNames))...)

    if 0 in sectors
        numsecs += 1
        replace!(sectors, 0 => numsecs)
    end
    
    sectorgroups = zeros(Int, numsecs)
    numsecgs = 1
    for iall in sp2
        i1 = iall[1]
        i2 = iall[2]
        if sectorgroups[sectors[i1]] == 0
            if sectorgroups[sectors[i2]] == 0
                sectorgroups[sectors[i1]] = 2 * numsecgs
                sectorgroups[sectors[i2]] = 2 * numsecgs + 1
                numsecgs += 1
            else
                sectorgroups[sectors[i1]] = sectorgroups[sectors[i2]] + 2 * iseven(sectorgroups[sectors[i2]]) - 1
            end
        else
            if sectorgroups[sectors[i2]] == 0
                sectorgroups[sectors[i2]] = sectorgroups[sectors[i1]] + 2 * iseven(sectorgroups[sectors[i1]]) - 1
            else
                sg1 = sectorgroups[sectors[i1]] ÷ 2
                sg2 = sectorgroups[sectors[i2]] ÷ 2
                if (sg1) != (sg2)
                    if sg1 < sg2
                        sectorgroups = map(
                            x -> begin
                                sgx = x ÷ 2
                                if sgx > sg2
                                    x - 2
                                elseif sgx == sg2
                                    2 * sg1 + x % 2
                                else
                                    x
                                end
                            end,
                            sectorgroups
                        )
                        # replace!(sectorgroups, 2 * sg2 + 0 => 2 * sg1 + 0)
                        # replace!(sectorgroups, 2 * sg2 + 1 => 2 * sg1 + 1)
                        numsecgs -= 1
                    else
                        sectorgroups = map(
                            x -> begin
                                sgx = x ÷ 2
                                if sgx > sg1
                                    x - 2
                                elseif sgx == sg1
                                    2 * sg2 + x % 2
                                else
                                    x
                                end
                            end,
                            sectorgroups
                        )
                        # replace!(sectorgroups, 2 * sg1 + 0 => 2 * sg2 + 0)
                        # replace!(sectorgroups, 2 * sg1 + 1 => 2 * sg2 + 1)
                        numsecgs -= 1
                    end
                end
            end
        end
    end
    sectorord = sortperm(sectorgroups)
    newsector = sort(eachindex(1:numsecs); by=i->sectorord[i])
    truesectors = map(j -> newsector[j], sectors)

    Psector = reduce(hcat, map(x -> map(y -> (x == y)+0, 1:D), sortperm(truesectors)))
    permBA = replacetags(eigenBA, "Eigen" => "Perm")
    permAB = replacetags(permBA, "BA", "AB")

    groupedSect = accumulate((x,y) -> x[1] == y ? (x[1], x[2]+1) : (x[1]+1, 1), sort(sectorgroups); init=(1,0))
    groupedSect = map(x -> (x[1] ÷ 2, x[1] % 2, x[2]), groupedSect)
    sectoredInd = accumulate((x,y) -> x[1] == y ? (x[1], x[2]+1) : (x[1]+1, 1), sort(truesectors); init=(1,0))
    sectoredInd = map(x -> begin
        y = groupedSect[x[1]]
        (y[1], y[2], y[3], x[2])
    end, sectoredInd)

    return ITensor(Psector, eigenBA, permBA), permBA, permAB, sectoredInd
end

function logSDG(target::String)
    f = h5open("./snapshots/$(target)/mps.h5")
    stA = read(f, "A", ITensor)
    stB = read(f, "B", ITensor)
    bwAB = read(f, "AB", ITensor)
    bwBA = read(f, "BA", ITensor)

    A, B, sA, sB, AB, BA = divideBondWeights(stA, stB, bwAB, bwBA)

    for is in eachval(sA)
        A_sg, B_sg, lBA, lAB = svdGauge(A, B, is, sA, BA)
        A_sdg, B_sdg, eBA, eAB = svdDiagGauge(A_sg, B_sg, lBA, lAB, is, sA, sB)
        Psec, pBA, pAB, sI = getSectorPerm(B_sdg, is, sB, eBA)
        A_sdgs = Psec * A_sdg * replaceinds(Psec, [pBA, eBA], [pAB, eAB])
        B_sdgs = replaceinds(Psec, [pBA, eBA], [pAB, eAB]) * B_sdg * Psec

        snapshot(target, "As$(is)", A_sg, sA, lAB, lBA)
        snapshot(target, "Bs$(is)", B_sg, sB, lAB, lBA)
        snapshot(target, "Asd$(is)", A_sdgs, sA, pAB, pBA)
        snapshot(target, "Bsd$(is)", B_sdgs, sB, pAB, pBA)
        logTM(target, "sd$(is)", A_sdgs, B_sdgs, pBA, sI)
    end
end



# logSDG("2024-07-16/Hubbard/U=10.0/iTEBD/mpslen=2/D=64/seed=10/initΔτ=0.1")

# logSDG("2024-07-02/Hubbard/U=10.0/iTEBD/mpslen=2/D=16/seed=10/initΔτ=0.1")
# logSDG("2024-07-02/Hubbard/U=10.0/iTEBD/mpslen=2/D=16/seed=10/initΔτ=0.1/TI")
# logSDG("2024-07-02/Hubbard/U=10.0/iTEBD/mpslen=2/D=16/seed=10/initΔτ=0.1/Mirror")
# logSDG("2024-07-02/Hubbard/U=10.0/iTEBD/mpslen=2/D=16/seed=10/initΔτ=0.1/MirrorTI")
# logSDG("2024-07-03/Hubbard/U=10.0/iTEBD/mpslen=2/D=16/seed=20/initΔτ=0.1")
# logSDG("2024-07-03/Hubbard/U=10.0/iTEBD/mpslen=2/D=16/seed=20/initΔτ=0.1/TI")
# logSDG("2024-07-03/Hubbard/U=10.0/iTEBD/mpslen=2/D=16/seed=20/initΔτ=0.1/Mirror")
# logSDG("2024-07-03/Hubbard/U=10.0/iTEBD/mpslen=2/D=16/seed=20/initΔτ=0.1/MirrorTI")
# logSDG("2024-07-03/Hubbard/U=10.0/iTEBD/mpslen=2/D=16/seed=30/initΔτ=0.1")
# logSDG("2024-07-03/Hubbard/U=10.0/iTEBD/mpslen=2/D=16/seed=30/initΔτ=0.1/TI")
# logSDG("2024-07-03/Hubbard/U=10.0/iTEBD/mpslen=2/D=16/seed=30/initΔτ=0.1/Mirror")
# logSDG("2024-07-03/Hubbard/U=10.0/iTEBD/mpslen=2/D=16/seed=30/initΔτ=0.1/MirrorTI")

# logSDG("2024-07-02/Hubbard/U=0.0/iTEBD/mpslen=2/D=16/seed=10/initΔτ=0.1")
# logSDG("2024-07-02/Hubbard/U=0.0/iTEBD/mpslen=2/D=16/seed=10/initΔτ=0.1/TI")
# logSDG("2024-07-02/Hubbard/U=0.0/iTEBD/mpslen=2/D=16/seed=10/initΔτ=0.1/Mirror")
# logSDG("2024-07-02/Hubbard/U=0.0/iTEBD/mpslen=2/D=16/seed=10/initΔτ=0.1/MirrorTI")
# logSDG("2024-07-03/Hubbard/U=0.0/iTEBD/mpslen=2/D=16/seed=20/initΔτ=0.1")
# logSDG("2024-07-03/Hubbard/U=0.0/iTEBD/mpslen=2/D=16/seed=20/initΔτ=0.1/TI")
# logSDG("2024-07-03/Hubbard/U=0.0/iTEBD/mpslen=2/D=16/seed=20/initΔτ=0.1/Mirror")
# logSDG("2024-07-03/Hubbard/U=0.0/iTEBD/mpslen=2/D=16/seed=20/initΔτ=0.1/MirrorTI")
# logSDG("2024-07-03/Hubbard/U=0.0/iTEBD/mpslen=2/D=16/seed=30/initΔτ=0.1")
# logSDG("2024-07-03/Hubbard/U=0.0/iTEBD/mpslen=2/D=16/seed=30/initΔτ=0.1/TI")
# logSDG("2024-07-03/Hubbard/U=0.0/iTEBD/mpslen=2/D=16/seed=30/initΔτ=0.1/Mirror")
# logSDG("2024-07-03/Hubbard/U=0.0/iTEBD/mpslen=2/D=16/seed=30/initΔτ=0.1/MirrorTI")

# logSDG("2024-07-02/Hubbard/U=-10.0/iTEBD/mpslen=2/D=16/seed=10/initΔτ=0.1")
# logSDG("2024-07-02/Hubbard/U=-10.0/iTEBD/mpslen=2/D=16/seed=10/initΔτ=0.1/TI")
# logSDG("2024-07-02/Hubbard/U=-10.0/iTEBD/mpslen=2/D=16/seed=10/initΔτ=0.1/Mirror")
# logSDG("2024-07-02/Hubbard/U=-10.0/iTEBD/mpslen=2/D=16/seed=10/initΔτ=0.1/MirrorTI")
# logSDG("2024-07-03/Hubbard/U=-10.0/iTEBD/mpslen=2/D=16/seed=20/initΔτ=0.1")
# logSDG("2024-07-03/Hubbard/U=-10.0/iTEBD/mpslen=2/D=16/seed=20/initΔτ=0.1/TI")
# logSDG("2024-07-03/Hubbard/U=-10.0/iTEBD/mpslen=2/D=16/seed=20/initΔτ=0.1/Mirror")
# logSDG("2024-07-03/Hubbard/U=-10.0/iTEBD/mpslen=2/D=16/seed=20/initΔτ=0.1/MirrorTI")
# logSDG("2024-07-03/Hubbard/U=-10.0/iTEBD/mpslen=2/D=16/seed=30/initΔτ=0.1")
# logSDG("2024-07-03/Hubbard/U=-10.0/iTEBD/mpslen=2/D=16/seed=30/initΔτ=0.1/TI")
# logSDG("2024-07-03/Hubbard/U=-10.0/iTEBD/mpslen=2/D=16/seed=30/initΔτ=0.1/Mirror")
# logSDG("2024-07-03/Hubbard/U=-10.0/iTEBD/mpslen=2/D=16/seed=30/initΔτ=0.1/MirrorTI")

for U in [10.0, 0.0, -10.0]
    for seed in [10, 20, 30]
        for type in ["/NoSymm", "/TI", "/Mirror", "/MirrorTI", "/SII", "/PHS", "/PHESII"]
            target = "2024-07-24/Hubbard/U=$(U)/iTEBD/mpslen=2/D=16/seed=$(seed)/initΔτ=0.1$(type)"
            println(target)
            logSDG(target)
        end
    end
end

# run(`export TARGET=$(target)`)

