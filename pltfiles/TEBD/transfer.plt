reset
target = system("echo $TARGET")

tmfile = "./snapshots/".target."/Asd1_tm.dat"
bonddim = system("tail -n 1 ".tmfile." | awk -F ',' '{print $1}'") + 0

set term png size 1600,1600

set output "./plots/".target."/tm.png"
set title target
set isotropic
unset key
set cbrange [0:1]
set xrange [0.5:0.5+bonddim]
set yrange [0.5:0.5+bonddim]
plot tmfile u 1:2:3 sparse matrix=(bonddim,bonddim) origin=(1,1) with image
unset output
