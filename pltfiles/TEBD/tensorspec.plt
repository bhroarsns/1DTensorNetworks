target = "F_XXX/S=1_2/iTEBD/mpslen=2/D=16/seed=10"
titlestr = "F XXX, S=1/2, iTEBD, mpslen=2, D=16, seed=10"
specfile = "./snapshots/".target."/spec.dat"

totsize = system("wc -l ".specfile." | awk '{print $1}'")
array SA[totsize]
array NA[totsize]
array RA[totsize]
array IA[totsize]
array AA[totsize]
array GA[totsize]
stats specfile u (SA[1+$0]=$1, NA[1+$0]=sprintf("%d", $2), RA[1+$0]=$3, IA[1+$0]=$4, AA[1+$0]=$5, GA[1+$0]=$6) nooutput

set term png size 1600,1600
unset key
set title titlestr
set xrange [-2.5:2.5]
set xtics ("↓" -2, "↑" -1, "↑" 1, "↓" 2)
# set xrange [-3.5:3.5]
# set xtics ("-" -3, "0" -2, "+" -1, "+" 1, "0" 2, "-" 3)
# set xrange [-4.5:4.5]
# set xtics ("UpDn" -4, "Dn" -3, "Up" -2, "Emp" -1, "Emp" 1, "Up" 2, "Dn" 3, "UpDn" 4)

outfile = "./plots/".target."/tensorspec.png"
set output outfile
do for [i=1:totsize] {
    if (IA[i] < 0.0) {
        set label NA[i] at (SA[i]+0.25*(IA[i]/AA[i])),(RA[i]) right font ",10"
    } else if (IA[i] > 0.0) {
        set label NA[i] at (SA[i]+0.25*(IA[i]/AA[i])),(RA[i]) left font ",10"
    } else {
        set label NA[i] at (SA[i]),(RA[i]) center font ",10"
    }
}
plot specfile u ($1+0.25*($4/$5)):3 ps 3

unset output
unset label

outfile = "./plots/".target."/tensorspecabs.png"
set output outfile
set yrange [-0.1:]
do for [i=1:totsize] {
    if (IA[i] < 0.0) {
        set label NA[i] at (SA[i]+0.25*(IA[i]/AA[i])),(AA[i]) right font ",10"
    } else if (IA[i] > 0.0) {
        set label NA[i] at (SA[i]+0.25*(IA[i]/AA[i])),(AA[i]) left font ",10"
    } else {
        set label NA[i] at SA[i],(AA[i]) center font ",10"
    }
}
plot specfile u ($1+0.25*($4/$5)):5 ps 3
unset output
