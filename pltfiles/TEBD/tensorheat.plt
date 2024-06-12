target = system("echo $TARGET")
sitetype = 4

filename(gt,ig,it,is) = sprintf("./snapshots/%s/%c%s%d_%d.dat", target, 64+it, gt == 1 ? "s" : "sd", ig, is)

# 各行列のサイズを1×1とする
# 行列間の余白幅をmとする
# 各gaugeのパネルのサイズは(st+(st-1)*m)×(2+m)
# gaugeのラベルは(st+(st-1)*m)×1ぐらい?
# gauge一個分で(st+(st-1)*m)×(3+m)
# gaugeパネル間の余白幅はMP*m
# →全体のサイズは(st*st+(st*st-st+(st+1)*MP)*m)×(6+(2+(st+1)*MP)*m)
mult = 10
MP = 4
M = 16
m = 1.0 / M
xlen = sitetype * sitetype + (sitetype * sitetype - sitetype + (sitetype + 1) * MP) * m
ylen = 6 + (2 + (sitetype + 1) * MP) * m
xsize = mult * M * xlen
ysize = mult * M * ylen
smatx = 1.0 / xlen
smarx = 1.0 / xlen * m
smaty = 1.0 / ylen
smary = 1.0 / ylen * m
sgaux = sitetype * smatx + (sitetype - 1) * smarx
sgauy = 3 * smaty + smary

originx(gt,ig,it,is) = MP * smarx + (ig - 1) * (sgaux + MP * smarx) + (is - 1) * (smatx + smarx)
originy(gt,ig,it,is) = MP * smary + (2 - gt) * (sgauy + MP * smary) + (2 - it) * (smaty + smary)

set term png size xsize,ysize
set output "./plots/".target."/heat.png"

set multiplot

set origin originx(0,2,0,1),(sgauy+4*MP*smary+smary+2.5*smaty)
set size 2*sgaux,smaty*0.5
set margins 0,0,0,0
unset tics
unset key
set xrange [1:3]
set yrange [-3:-1]
set label target at 2,-2 center font ",24" noenhanced
plot '+'
unset label

do for [ig = 1:sitetype] {
    set origin originx(0,ig,0,1),(2*MP*smary+smary+2*smaty)
    set size sgaux,smaty*0.5
    set margins 0,0,0,0
    unset tics
    unset key
    set xrange [1:3]
    set yrange [-3:-1]
    set label sprintf("A%ds-B%dd gauge", ig, sitetype - ig + 1) at 2,-2 center font ",24"
    plot '+'
    unset label

    set origin originx(0,ig,0,1),(sgauy+3*MP*smary+smary+2*smaty)
    set size sgaux,smaty*0.5
    set margins 0,0,0,0
    unset tics
    unset key
    set xrange [1:3]
    set yrange [-3:-1]
    set label sprintf("A%ds gauge", ig) at 2,-2 center font ",24"
    plot '+'
    unset label
}

set border
do for [gt=1:2] {
    do for [ig=1:sitetype] {
        do for [it=1:2] {
            do for [is=1:sitetype] {
                unset tics
                unset key
                set origin originx(gt,ig,it,is),originy(gt,ig,it,is)
                set size smatx,smaty
                set margins 0,0,0,0
                file = filename(gt,ig,it,is)
                DAB = system("awk 'NR==1{print $1}' ".file) + 0.0
                DBA = system("awk 'NR==1{print $2}' ".file) + 0.0
                totsize = DAB * DBA
                set xrange [0.5:(DAB+0.5)]
                set yrange [0.5:(DBA+0.5)]
                array XC[totsize]
                array YC[totsize]
                array ZC[totsize]
                stats filename(gt,ig,it,is) every ::1 u (XC[$0+1] = $1, YC[$0+1] = $2, ZC[$0+1] = $3) nooutput
                unset cbtics
                set cbrange [0:1]
                unset colorbox
                plot ZC u (XC[$1]):(YC[$1]):(ZC[$1]) sparse matrix=(DAB,DBA) origin=(1,1) with image
                # plot '+'
            }
        }
    }
}

unset multiplot
unset output
