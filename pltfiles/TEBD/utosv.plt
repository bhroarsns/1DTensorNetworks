D = 16
seed = 10
files = system(sprintf("ls -d ./snapshots/2024-06-05/Hubbard/U=*/iTEBD/mpslen=2/D=%d/seed=%d/initΔτ=0.1/", D, seed))
files = split(files, "\n")
nfile = |files|

set xrange [-10.0:10.0]
set yrange [0.0:1.0]

set term png size 6400,3200
set output "plots/2024-06-05/Hubbard/test.png"
set multiplot
unset key

do for [is=1:4] {
    do for [it=0:1] {
        set size 0.25,0.5
        set origin (is-1)*0.25,it*0.5
        do for [id=1:D] {
            print(sprintf("state #%d, tensor #%c, singular value #%d", is, it+65, id))
            array Uarr[nfile]
            array Varr[nfile]
            do for [ifile=1:nfile] {
                file = sprintf("%s%c%d.dat", files[ifile], it+65, is)
                Uarr[ifile] = system("echo ".file." | grep -oE \"U=[\-]*[0-9\.]+\" | sed \"s/U=//g\"")
                Varr[ifile] = system(sprintf("awk -F ', ' 'END{print $%d}' %s", id, file))
            }
            plot Uarr u (Uarr[$1]+0):(Varr[$1]+0) w l smooth unique
        }
    }
}
unset multiplot
unset output


# files = system(sprintf("ls ./snapshots/Hubbard/U=*/iTEBD/mpslen=2/D=%d/seed=10/singular.dat", D))
# files = split(files, "\n")
# nfile = |files|

# set xrange [-10.0:10.0]
# set yrange [0.0:1.0]

# array Uarr[D*nfile]
# array Varr[D*nfile]
# do for [id=1:D] {
#     do for [ifile=1:nfile] {
#         file = files[ifile]
#         Uarr[nfile*(id-1)+ifile] = system("echo ".file." | grep -oE \"U=[\-]*[0-9\.]+\" | sed \"s/U=//g\"")
#         Varr[nfile*(id-1)+ifile] = system(sprintf("awk -F ', ' '$1 == \"%d\" && $2 == \"%d\" { print $3 }' %s", 1, id, file))
#     }
# }
# plot for [id=1:D] [(id-1)*nfile+1:id*nfile] (Uarr[$x]+0),(Varr[$x]+0)
# pause -1

