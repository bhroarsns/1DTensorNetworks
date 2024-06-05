set term png
target = "2024-06-05/Hubbard/U=0.0/iTEBD/mpslen=2/D=16/seed=10/initΔτ=0.1"

set xrange [0:0.1]
set output "plots/".target."/order.png"
# plot "results/".target."/energy.dat" u (1/$2):(abs($7 - $8)) w l title "|<nA>-<nB>|" lw 4, "results/".target."/energy.dat" u (1/$2):(abs($10 - $13)) w l title "|<nA↑>-<nA↓>|" lw 4
plot "results/".target."/energy.dat" u (1/$2):(abs($9 - $10)) w l title "|<nA>-<nB>|" lw 4, "results/".target."/energy.dat" u (1/$2):(abs($12 - $15)) w l title "|<nA↑>-<nA↓>|" lw 4
unset output