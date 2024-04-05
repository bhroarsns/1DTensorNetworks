set term tikz standalone size 8in,6in
set output "./tex/F_XXX/S=1_2/iTEBD/mpslen=2/D=16/seed=10/energy.tex"

set xlabel 'TEBD steps'
set ylabel 'Energy Density $E_{GS}$'
set xrange [0:500]

plot "results/F_XXX/S=1_2/iTEBD/mpslen=2/D=16/seed=10/energy.dat" u 1:3 w l title "TEBD result", -0.25 w l title "$-0.25$"

set xlabel 'TEBD steps'
set ylabel 'Relative Error of Energy Density $\Delta E_{GS}$'
set logscale y
set xrange [0:500]

plot "results/F_XXX/S=1_2/iTEBD/mpslen=2/D=16/seed=10/energy.dat" u 1:(abs($3+0.25)/0.25) w l title "TEBD result"

set xlabel 'Inverse of TEBD steps'
set ylabel 'Relative Error of Energy Density $\Delta E_{GS}$'
unset logscale
set xrange [0:1.0]
plot "results/F_XXX/S=1_2/iTEBD/mpslen=2/D=16/seed=10/energy.dat" u (1/$1):(abs($3+0.25)/0.25) w l title "TEBD result"

set xrange [0:0.1]
plot "results/F_XXX/S=1_2/iTEBD/mpslen=2/D=16/seed=10/energy.dat" u (1/$1):(abs($3+0.25)/0.25) w l title "TEBD result"

unset output
system("/Library/TeX/texbin/pdflatex -output-directory ./plots/F_XXX/S=1_2/iTEBD/mpslen=2/D=16/seed=10 ./tex/F_XXX/S=1_2/iTEBD/mpslen=2/D=16/seed=10/energy.tex;")

system("echo latex")