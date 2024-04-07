set term tikz standalone size 8in,6in
set output "./tex/Hubbard/U=0.0/iTEBD/mpslen=2/D=32/seed=10/spectrum.tex"
file = "snapshots/Hubbard/U=0.0/iTEBD/mpslen=2/D=32/seed=10/spectrum/AB_left.dat"
complex(x,y) = x*{1,0}+y*{0,1}

set xlabel 'TEBD steps'
set ylabel 'Transfer Matrix Eigenvalues'
set xrange [0:500]

set yrange [0:2]
plot for [i=1:32] file u 0:(abs(complex(column(2*i-1),column(2*i)))) w l title sprintf("$|\\lambda_{%d}|$", i)

set yrange [0:1]
plot for [i=1:31] file u 0:(abs(complex(column(2*i+1),column(2*i+2)))/abs(complex(column(1),column(2)))) w l title sprintf("$|\\lambda_{%d}|/|\\lambda_1|$", i+1)

unset output
system("/Library/TeX/texbin/pdflatex -output-directory ./plots/Hubbard/U=0.0/iTEBD/mpslen=2/D=32/seed=10 ./tex/Hubbard/U=0.0/iTEBD/mpslen=2/D=32/seed=10/spectrum.tex;")

system("echo latex")