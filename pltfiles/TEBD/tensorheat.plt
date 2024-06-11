target = "2024-06-11/Hubbard/U=0.0/iTEBD/mpslen=2/D=32/seed=10/initΔτ=0.1"
sitetype = 4

sdir = "./snapshotdir/".target

# 各行列のサイズを1×1とする
# 行列間の余白幅をmとする
# 各gaugeのパネルのサイズは(st+(st-1)*m)×(2+m)
# gaugeのラベルは(st+(st-1)*m)×1ぐらい?
# gauge一個分で(st+(st-1)*m)×(3+m)
# gaugeパネル間の余白幅は2m
# →全体のサイズは(st*st+(st*st+st+2)*m)×(6+8m)
mult = 10
M = 16
m = 1.0 / M
xlen = sitetype * sitetype + (sitetype * sitetype + sitetype + 2) * m
ylen = (6 + 8 * m)
xsize = mult * M * xlen
ysize = mult * M * ylen
smatx = 1.0 / xlen
smarx = 1.0 / xlen * m
smaty = 1.0 / ylen
smary = 1.0 / ylen * m
sgaux = sitetype * smatx + (sitetype - 1) * smarx
sgauy = 3 * smaty + smary

originx(gt,ig,it,is) = smarx + (ig - 1) * (sgaux + 2 * smarx) + (is - 1) * (smatx + smarx)
originy(gt,ig,it,is) = 1 - 2 * smary - (gt - 1) * (sgauy + 2 * smary) - (it - 1) * (smaty + smary)

# set term png size xsize,ysize
# set output "heat.png"

set multiplot

do for [gt=1:2] {
    do for [ig=1:sitetype] {
        do for [it=1:2] {
            do for [is=1:sitetype] {
                unset tics
                unset key
                print(originy(gt,ig,it,is))
                set origin originx(gt,ig,it,is),originy(gt,ig,it,is)
                set size square smatx,smaty
                set xrange [-2:-1]
                set yrange [1:2]
                plot '+'
            }
        }
    }
}

unset multiplot
unset output

pause -1