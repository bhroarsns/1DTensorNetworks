reset
target = system("echo $TARGET")
outfile = "./tex/".target."/spectrum.tex"
sfile = "snapshots/".target."/Spec/FUN/AB_left_sym.dat"
afile = "snapshots/".target."/Spec/FUN/AB_left_asym.dat"
tfile = "snapshots/".target."/Spec/FUN/AB_left_tot.dat"

set term tikz standalone size 8in,6in
set output outfile
# set title target
complex(x,y) = x*{1,0}+y*{0,1}

set xlabel 'TEBD steps'
set ylabel 'Transfer Matrix Eigenvalues'
set xrange [0:1000]
unset key

set yrange [0:2]
plot for [i=1:16] sfile u 0:(abs(complex(column(2*i),column(2*i+1)))) w l lc 1, for [i=1:16] afile u 0:(abs(complex(column(2*i),column(2*i+1)))) w l lc 2

# set yrange [0:1]
# plot for [i=1:31] file u 0:(abs(complex(column(2*i+1),column(2*i+2)))/abs(complex(column(1),column(2)))) w l

unset output
system("/Library/TeX/texbin/pdflatex -output-directory ./plots/".target." ".outfile.";")

system("echo latex")