target = "Hubbard/U=0.0/iTEBD/mpslen=2/D=16/seed=12"
sitetype = 4
titlestr = system("echo ".target." | sed 's|\/|, |g' | sed 's|1_2|1\/2|' | sed 's|_|\ |g'")
specfile = "./snapshots/".target."/singular.dat"

totsize = system("wc -l ".specfile." | awk '{print $1}'")
array SA[totsize]
array NA[totsize]
array RA[totsize]
stats specfile u (SA[1+$0]=$1, NA[1+$0]=sprintf("%d", $2), RA[1+$0]=$3) nooutput

set term tikz standalone size 6in,6in
# set term png size 1600,1600
unset key
set title titlestr
set grid
if (sitetype == 2) {
    set xrange [-2.5:2.5]
    set xtics ("↓" -2, "↑" -1, "↑" 1, "↓" 2)
} else if (sitetype == 3) {
    set xrange [-3.5:3.5]
    set xtics ("-" -3, "0" -2, "+" -1, "+" 1, "0" 2, "-" 3)
} else if (sitetype == 4) {
    set xrange [-4.5:4.5]
    set xtics ("UpDn" -4, "Dn" -3, "Up" -2, "Emp" -1, "Emp" 1, "Up" 2, "Dn" 3, "UpDn" 4)
}

outfile = "./tex/".target."/tensorsing.tex"
set output outfile
do for [i=1:totsize] {
    if (NA[i] % 2 == 1) {
        set label NA[i] at SA[i]-0.15-0.35*(NA[i]/16.0),RA[i]+(0.5-((NA[i]-1)/2)%2)*0.01 left font ",8"
    } else {
        set label NA[i] at SA[i]+0.15+0.35*(NA[i]/16.0),RA[i]+(0.5-((NA[i]-1)/2)%2)*0.01 right font ",8"
    }
}
plot specfile u 1:3 ps 2 pt 6 lc 1, specfile u 1:3 ps 2 pt 1 lc 1

unset output
system("/Library/TeX/texbin/pdflatex -output-directory ./plots/".target." ".outfile.";")

# unset label
# unset title

# abfile = "./snapshots/".target."/AB.dat"
# bafile = "./snapshots/".target."/BA.dat"

# outfile = "./plots/".target."/tensorheat.png"
# totheight = 6400*sitetype
# set term png size 6400,4800
# set output outfile
# set margin 10, 10, 10, 10
# set multiplot
# unset key

# set size 1,0.1
# set origin 0,0.9
# set label titlestr font ",96" at 0,0 center
# unset border
# unset tics
# unset grid
# set xrange [-0.5:0.5]
# set yrange [-0.5:0.5]
# plot '+'
# unset label

# if (sitetype * (sitetype + 1) / 9.0 < 0.45) {
#     minia = 1.0 / sitetype * (sitetype + 1)
# } else {
#     minia = 0.45 / 9.0
# }
# miniw = minia * (sitetype + 1)
# minih = minia * 4.5

# do for [j=0:2*sitetype] {
#     add = (j==0 ? "" : ((j > sitetype) ? sprintf("%dsd",j-sitetype) : sprintf("%ds",j)))
#     afile = "./snapshots/".target."/A".add.".dat"
#     bfile = "./snapshots/".target."/B".add.".dat"

#     da = system("awk -F \", \" 'NR==1{print $1}' ".afile)
#     Dab = system("awk -F \", \" 'NR==1{print $2}' ".afile)
#     Dba = system("awk -F \", \" 'NR==1{print $3}' ".afile)
#     db = system("awk -F \", \" 'NR==1{print $1}' ".bfile)

#     sizea = da * Dba * Dab
#     sizeb = db * Dab * Dba
#     wid = da * (1 + Dba) - 1
#     hei = Dab
#     cx = 0
#     cy = hei + 1

#     sx(is, iab, iba) = iba + (is - 1) * (1 + Dba)
#     sy(is, iab, iba) = iab

#     offset0 = 1
#     offset1 = offset0 + sizea
#     # offset2 = offset1 + sizeb
#     # offset3 = offset2 + da*Dba
#     # offset4 = offset3 + Dab
#     # offset5 = offset4 + Dab
#     # totsize = sizea + sizeb + da*Dba + 2*Dab
#     totsize = sizea + sizeb

#     array XC[totsize]
#     array YC[totsize]
#     array ZC[totsize]

#     # array BBA[Dba]
#     # array BAB[Dab]
#     # stats bafile u (BBA[$1] = $2) nooutput
#     # stats abfile u (BAB[$1] = $2) nooutput

#     stats afile every ::1 u (XC[$0+offset0] = cx + sx($1, $2, $3), YC[$0+offset0] = cy + sy($1, $2, $3), ZC[$0+offset0] = $6) nooutput
#     stats bfile every ::1 u (XC[$0+offset1] = cx + sx($1, $2, $3), YC[$0+offset1] = cy - sy($1, $2, $3), ZC[$0+offset1] = $6) nooutput

#     # do for [i=1:da] {
#     #     stats bafile u (XC[$0+offset2+(i-1)*Dba] = cx + sx(i, 0, $1), YC[$0+offset2+(i-1)*Dba] = cy, ZC[$0+offset2+(i-1)*Dba] = $2) nooutput
#     # }
#     # stats abfile u (XC[$0+offset3] = cx, YC[$0+offset3] = cy + sy(1, $1, 0), ZC[$0+offset3] = $2) nooutput
#     # stats abfile u (XC[$0+offset4] = cx, YC[$0+offset4] = cy - sy(1, $1, 0), ZC[$0+offset4] = $2) nooutput

#     set isotropic

#     set xrange [0.5:0.5+wid]
#     set yrange [0.5:0.5+2*hei+1]
#     set margins 10, 10, 0,0
#     set ytics ("A" cy+(1+Dab)*0.5, "B" cy-(1+Dab)*0.5) font ",64" offset -2,0
#     if (sitetype == 2) {
#         set xtics ("↑" cx+(1+Dba)*0.5, "↓" cx+(1+Dba)*1.5) font ",64" offset 0,-6
#     } else if (sitetype == 3) {
#         set xtics ("+" cx+(1+Dba)*0.5, "0" cx+(1+Dba)*1.5, "-" cx+(1+Dba)*2.5) font ",64" offset 0,-6
#     } else if (sitetype == 4) {
#         set xtics ("Emp" cx+(1+Dba)*0.5, "Up" cx+(1+Dba)*1.5, "Dn" cx+(1+Dba)*2.5, "UpDn" cx+(1+Dba)*3.5) font ",64" offset 0,-6
#     }

#     if (j==0) {
#         set title "Canonical gauge" font ",64"
#         set size 0.5,0.45
#         set origin 0.25,0.45
#     } else {
#         if (j > sitetype) {
#             set title sprintf("A%ds-B%dd gauge", j - sitetype, 2*sitetype-j+1)
#         } else {
#             set title sprintf("A%ds gauge", j) font ",64"
#         }
#         xwidth = 1.0 / sitetype
#         minix = ((j-1)%sitetype)*xwidth + (xwidth - miniw) / 2.0
#         miniy = 0.225 - 0.225 * ((j-1) / sitetype)
#         set size miniw,minih
#         set origin minix,miniy
#     }
#     set cbtics
#     set cbrange [0:1]
#     plot ZC u (XC[$1]):(YC[$1]):(ZC[$1]) sparse matrix=(wid,2*hei+1) origin=(1,1) with image
# }

# unset multiplot
# # pause -1

# unset output
