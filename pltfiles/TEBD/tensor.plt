target = "Hubbard/U=10.0/iTEBD/mpslen=2/D=16/seed=10"
titlestr = "Hubbard, U=10.0, iTEBD, mpslen=2, D=16, seed=10"
abfile = "./snapshots/".target."/AB.dat"
bafile = "./snapshots/".target."/BA.dat"

do for [j=0:4] {
    add = (j==0 ? "" : sprintf("%d",j))
    outfile = "./plots/".target."/tensor".add.".png"
    afile = "./snapshots/".target."/A".add.".dat"
    bfile = "./snapshots/".target."/B".add.".dat"

    da = system("awk -F \", \" 'NR==1{print $1}' ".afile)
    Dab = system("awk -F \", \" 'NR==1{print $2}' ".afile)
    Dba = system("awk -F \", \" 'NR==1{print $3}' ".afile)
    db = system("awk -F \", \" 'NR==1{print $1}' ".bfile)

    sizea = da * Dba * Dab
    sizeb = db * Dab * Dba
    wid = da * (1 + Dba)
    hei = 1 + Dab
    cx = 1
    cy = hei + 1

    sx(is, iab, iba) = 1 + iba + (is - 1) * (1 + Dba)
    sy(is, iab, iba) = 1 + iab

    offset0 = 1
    offset1 = offset0 + sizea
    offset2 = offset1 + sizeb
    offset3 = offset2 + da*Dba
    offset4 = offset3 + Dab
    offset5 = offset4 + Dab
    totsize = sizea + sizeb + da*Dba + 2*Dab

    array XC[totsize]
    array YC[totsize]
    array ZC[totsize]

    # array BBA[Dba]
    # array BAB[Dab]
    # stats bafile u (BBA[$1] = $2) nooutput
    # stats abfile u (BAB[$1] = $2) nooutput

    stats afile every ::1 u (XC[$0+offset0] = cx + sx($1, $2, $3), YC[$0+offset0] = cy + sy($1, $2, $3), ZC[$0+offset0] = $6) nooutput
    stats bfile every ::1 u (XC[$0+offset1] = cx + sx($1, $2, $3), YC[$0+offset1] = cy - sy($1, $2, $3), ZC[$0+offset1] = $6) nooutput

    do for [i=1:da] {
        stats bafile u (XC[$0+offset2+(i-1)*Dba] = cx + sx(i, 0, $1), YC[$0+offset2+(i-1)*Dba] = cy, ZC[$0+offset2+(i-1)*Dba] = $2) nooutput
    }
    stats abfile u (XC[$0+offset3] = cx, YC[$0+offset3] = cy + sy(1, $1, 0), ZC[$0+offset3] = $2) nooutput
    stats abfile u (XC[$0+offset4] = cx, YC[$0+offset4] = cy - sy(1, $1, 0), ZC[$0+offset4] = $2) nooutput

    set term png size 3200,2000
    set output outfile
    set isotropic
    unset key

    set title titlestr font ",56"
    set xrange [0.5:0.5+wid+1]
    set margins 10, 10, 10, 10
    set ytics ("A" cy+1+(1+Dab)*0.5, "B" cy-1-(1+Dab)*0.5) font ",96" offset -2,0
    set xtics ("Emp" cx+1+(1+Dba)*0.5, "Up" cx+1+(1+Dba)*1.5, "Dn" cx+1+(1+Dba)*2.5, "UpDn" cx+1+(1+Dba)*3.5) font ",96" offset 0,-6
    # set xtics ("+" cx+1+(1+Dba)*0.5, "0" cx+1+(1+Dba)*1.5, "-" cx+1+(1+Dba)*2.5) font ",96" offset 0,-6
    # set xtics ("↑" cx+1+(1+Dba)*0.5, "↓" cx+1+(1+Dba)*1.5) font ",96" offset 0,-6
    set yrange [0.5:0.5+2*hei+1]

    plot ZC u (XC[$1]):(YC[$1]):(ZC[$1]) sparse matrix=(wid+1,2*hei+1) origin=(1,1) with image

    unset output
}