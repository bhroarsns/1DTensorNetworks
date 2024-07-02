using ITensors
using HDF5

function snapshot(targetDir::String, name::String, tensor::ITensor, s::Index, bondAB::Index, bondBA::Index)
    for is in eachval(s)
        open("$(targetDir)/$(name)_$(is).dat", "w") do io
            println(io, dim(bondAB), " ", dim(bondBA))
            for ib in eachindval(bondAB, bondBA)
                entry = tensor[s=>is, ib...]
                println(io, join(getfield.(ib, 2), ", "), ", ", abs(entry), ", ", angle(entry), ", ", real(entry), ", ", imag(entry))
            end
        end
    end
end

target = "./snapshots/2024-07-02/Hubbard/U=0.0/iTEBD/mpslen=2/D=16/seed=10/initΔτ=0.1/MirrorTI"

f = h5open("$(target)/mps.h5")
stA = read(f, "A", ITensor)
stB = read(f, "B", ITensor)
bwAB = read(f, "AB", ITensor)
bwBA = read(f, "BA", ITensor)

A_BA, sA, A_AB = inds(stA)
B_AB, sB, B_BA = inds(stB)

AB = removetags(sim(A_AB), "A")
BA = removetags(sim(A_BA), "A")

A = replaceinds(sqrt.(bwBA) * stA * sqrt.(bwAB), [B_AB, B_BA], [AB, BA])
B = replaceinds(sqrt.(bwAB) * stB * sqrt.(bwBA), [A_AB, A_BA], [AB, BA])

for is in eachval(sA)
    matA = A * onehot(sA => is)
    U, S, V, _, linku, linkv = svd(matA, [BA])
    A_sg = dag(U) * A * dag(V)
    B_sg = V * B * U

    matB = B_sg * onehot(sB => dim(sA) - is + 1)
    D, P, _, e2, e = eigen(matB, [linkv], [linku])
    Pinv = ITensor(inv(matrix(P)), e2, linkv)

    A_sdg = replaceinds(Pinv, [e, linku], [e2, linkv]) * A_sg * replaceinds(P, [e2, linkv], [e, linku])
    B_sdg = Pinv * B_sg * P

    snapshot(target, "As$(is)", A_sg, sA, linkv, linku)
    snapshot(target, "Bs$(is)", B_sg, sB, linkv, linku)
    snapshot(target, "Asd$(is)", A_sdg, sA, e2, e)
    snapshot(target, "Bsd$(is)", B_sdg, sB, e2, e)
end