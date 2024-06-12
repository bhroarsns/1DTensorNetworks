set term png
target = system("echo $TARGET")

set xrange [0:0.5]
set logscale y
set key right bottom
set output "plots/".target."/order.png"
plot "results/".target."/energy.dat" u (1/$2):(abs($7 - $8)) w l title "|<nA>-<nB>|" lw 4, "results/".target."/energy.dat" u (1/$2):(abs($10 - $13)) w l title "|<nA↑>-<nA↓>|" lw 4, "results/".target."/energy.dat" u (1/$2):(abs($9 - $12)) w l title "|<n↑>-<n↓>|" lw 4
# plot "results/".target."/energy.dat" u (1/$2):(abs($9 - $10)) w l title "|<nA>-<nB>|" lw 4, "results/".target."/energy.dat" u (1/$2):(abs($12 - $15)) w l title "|<nA↑>-<nA↓>|" lw 4, "results/".target."/energy.dat" u (1/$2):(abs($11 - $14)) w l title "|<n↑>-<n↓>|" lw 4
unset output