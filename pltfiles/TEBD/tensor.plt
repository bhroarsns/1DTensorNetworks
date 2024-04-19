target = "Hubbard/U=-1.0/iTEBD/mpslen=2/D=16/seed=10"
outfile = "./plots/".target."/tensor.png"
afile = "./snapshots/".target."/Step/500/FUN/A"
bfile = "./snapshots/".target."/Step/500/FUN/B"
abfile = "./snapshots/".target."/Step/500/FUN/AB"
bafile = "./snapshots/".target."/Step/500/FUN/BA"

da = system("awk -F \", \" 'NR==2{print $1}' ".afile)
Dba = system("awk -F \", \" 'NR==2{print $2}' ".afile)
Dab = system("awk -F \", \" 'NR==2{print $3}' ".afile)
db = system("awk -F \", \" 'NR==2{print $1}' ".bfile)

set term png size 6400,2000
set output outfile
set isotropic
unset key

sizea = da * Dba * Dab
sizeb = db * Dab * Dba
wid = da * (1 + Dba)
hei = 1 + Dab
cx = wid + 1
cy = hei + 1

set xrange [0.5:0.5+2*wid+1]
set yrange [0.5:0.5+2*hei+1]

sx(is, iba, iab) = 1 + iba + (is - 1) * (1 + Dba)
sy(is, iba, iab) = 1 + iab

offset0 = 1
offset1 = offset0 + sizea
offset2 = offset1 + sizeb
offset3 = offset2 + sizea
offset4 = offset3 + sizeb
offset5 = offset4 + Dba
offset6 = offset5 + Dba
offset7 = offset6 + Dba
offset8 = offset7 + Dba
offset9 = offset8 + Dab
offset10 = offset9 + Dba
offset11 = offset10 + Dba
offset12 = offset11 + Dba
offset13 = offset12 + Dba
totsize = 2*sizea + 2*sizeb + 8*Dba + 2*Dab

array XC[totsize]
array YC[totsize]
array ZC[totsize]

array BBA[Dba]
array BAB[Dab]

stats bafile u (BBA[$1] = $2) nooutput
stats abfile u (BAB[$1] = $2) nooutput

stats afile every ::1 u (XC[$0+offset0] = cx + sx($1, $2, $3), YC[$0+offset0] = cy + sy($1, $2, $3), ZC[$0+offset0] = $6 * sqrt(BBA[$2]) * sqrt(BAB[$3])) nooutput
stats bfile every ::1 u (XC[$0+offset1] = cx + sx($1, $3, $2), YC[$0+offset1] = cy - sy($1, $3, $2), ZC[$0+offset1] = $6 * sqrt(BBA[$3]) * sqrt(BAB[$2])) nooutput
stats afile every ::1 u (XC[$0+offset2] = cx - sx($1, $2, $3), YC[$0+offset2] = cy - sy($1, $2, $3), ZC[$0+offset2] = $6 * sqrt(BBA[$2]) * sqrt(BAB[$3])) nooutput
stats bfile every ::1 u (XC[$0+offset3] = cx - sx($1, $3, $2), YC[$0+offset3] = cy + sy($1, $3, $2), ZC[$0+offset3] = $6 * sqrt(BBA[$3]) * sqrt(BAB[$2])) nooutput

stats bafile u (XC[$0+offset4]  = cx + sx(1, $1, 0), YC[$0+offset4]  = cy,                ZC[$0+offset4]  = $2) nooutput
stats bafile u (XC[$0+offset5]  = cx + sx(2, $1, 0), YC[$0+offset5]  = cy,                ZC[$0+offset5]  = $2) nooutput
stats bafile u (XC[$0+offset6]  = cx + sx(3, $1, 0), YC[$0+offset6]  = cy,                ZC[$0+offset6]  = $2) nooutput
stats bafile u (XC[$0+offset7]  = cx + sx(4, $1, 0), YC[$0+offset7]  = cy,                ZC[$0+offset7]  = $2) nooutput
stats abfile u (XC[$0+offset8]  = cx,                YC[$0+offset8]  = cy + sy(1, 0, $1), ZC[$0+offset8]  = $2) nooutput
stats bafile u (XC[$0+offset9]  = cx - sx(1, $1, 0), YC[$0+offset9]  = cy,                ZC[$0+offset9]  = $2) nooutput
stats bafile u (XC[$0+offset10] = cx - sx(2, $1, 0), YC[$0+offset10] = cy,                ZC[$0+offset10] = $2) nooutput
stats bafile u (XC[$0+offset11] = cx - sx(3, $1, 0), YC[$0+offset11] = cy,                ZC[$0+offset11] = $2) nooutput
stats bafile u (XC[$0+offset12] = cx - sx(4, $1, 0), YC[$0+offset12] = cy,                ZC[$0+offset12] = $2) nooutput
stats abfile u (XC[$0+offset13] = cx,                YC[$0+offset13] = cy - sy(1, 0, $1), ZC[$0+offset13] = $2) nooutput

set cbrange [0:1]

plot ZC u (XC[$1]):(YC[$1]):(ZC[$1]) sparse matrix=(2*wid+1,2*hei+1) origin=(1,1) with image

unset output