target = "Hubbard/U=10.0/iTEBD/mpslen=2/D=16/seed=10"
datfile = "results/".target."/energy.dat"
outfile = "tex/".target."/energy.tex"
exac = -0.267155

set term tikz standalone size 8in,6in
set output outfile
# set title target

set xlabel 'TEBD steps'
set ylabel 'Energy Density $E_{GS}$'
set xrange [0:500]

plot datfile u 1:2 w l title "TEBD result", exac w l title "Exact"

set xlabel 'TEBD steps'
set ylabel 'Relative Error of Energy Density $\Delta E_{GS}$'
set logscale y
set xrange [0:500]

plot datfile u 1:(abs($2-exac)/abs(exac)) w l title "TEBD result"

set xlabel 'Temperature $T$'
set ylabel 'Relative Error of Energy Density $\Delta E_{GS}$'
unset logscale
set xrange [0:1.0]
set yrange [0:]
plot datfile u (1/($1*0.1)):(abs($2-exac)/abs(exac)) w l title "TEBD result"

set xrange [0:0.1]
plot datfile u (1/($1*0.1)):(abs($2-exac)/abs(exac)) w l title "TEBD result"

unset output
system("/Library/TeX/texbin/pdflatex -output-directory ./plots/".target." ".outfile.";")

system("echo latex")