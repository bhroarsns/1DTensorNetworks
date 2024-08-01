function printCorrLabel(target)
    files = split(chop(read(`ls snapshots/$target/Corr`, String)), '\n')
    vals = map(file -> parse.(Float64, split(replace(read(`tail -n 1 snapshots/$target/Corr/$file`, String), '\n' => "", ' ' => ""), ',')), files)
    sp = sortperm(vals; by=x->x[2])
    group = accumulate((x,y) -> vals[y][2] < sqrt(eps()) ? x : isapprox(vals[x][2], vals[y][2]) ? x : y, sp; init=sp[1])
    ug = unique(group)

    maxval = vals[ug[end]][2]
    shift = accumulate((x, y) -> max(vals[y][2], x + maxval*0.01), ug; init=-0.01)

    open("snapshots/$target/CorrLabel.txt", "w") do io
        for (ig, g) in Iterators.enumerate(ug)
            print(io, vals[g][1], ", ", vals[g][2], ", ", shift[ig], ", ")
            for ind in findall(isequal(g), group)
                print(io, replace(files[sp[ind]], ".dat" => ""), " ")
            end
            println(io, "")
        end
    end
end

for U in [10.0, 0.0, -10.0]
    for seed in [10, 20, 30]
        for type in ["/NoSymm", "/TI", "/Mirror", "/MirrorTI", "/SII", "/PHS", "/PHESII"]
            target = "2024-07-24/Hubbard/U=$(U)/iTEBD/mpslen=2/D=16/seed=$(seed)/initΔτ=0.1$(type)"
            println(target)
            printCorrLabel(target)
        end
    end
end

# printCorrLabel("2024-07-16/Hubbard/U=10.0/iTEBD/mpslen=2/D=64/seed=10/initΔτ=0.1")