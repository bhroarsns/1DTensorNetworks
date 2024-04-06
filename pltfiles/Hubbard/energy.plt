set term tikz standalone size 8in,6in
set output "./tex/Hubbard/U=0.0/iTEBD/mpslen=2/D=16/seed=10/energy.tex"

set xlabel 'TEBD steps'
set ylabel 'Energy Density $E_{GS}$'
set xrange [0:500]

plot "results/Hubbard/U=0.0/iTEBD/mpslen=2/D=16/seed=10/energy.dat" u 1:3 w l title "TEBD result", -1.27324 w l title "$-1.27324$"

set xlabel 'TEBD steps'
set ylabel 'Relative Error of Energy Density $\Delta E_{GS}$'
set logscale y
set xrange [0:500]

plot "results/Hubbard/U=0.0/iTEBD/mpslen=2/D=16/seed=10/energy.dat" u 1:(abs($3+1.27324)/1.27324) w l title "TEBD result"

set xlabel 'Inverse of TEBD steps'
set ylabel 'Relative Error of Energy Density $\Delta E_{GS}$'
unset logscale
set xrange [0:1.0]
plot "results/Hubbard/U=0.0/iTEBD/mpslen=2/D=16/seed=10/energy.dat" u (1/$1):(abs($3+1.27324)/1.27324) w l title "TEBD result"

set xrange [0:0.1]
plot "results/Hubbard/U=0.0/iTEBD/mpslen=2/D=16/seed=10/energy.dat" u (1/$1):(abs($3+1.27324)/1.27324) w l title "TEBD result"

unset output
system("/Library/TeX/texbin/pdflatex -output-directory ./plots/Hubbard/U=0.0/iTEBD/mpslen=2/D=16/seed=10 ./tex/Hubbard/U=0.0/iTEBD/mpslen=2/D=16/seed=10/energy.tex;")

system("echo latex")