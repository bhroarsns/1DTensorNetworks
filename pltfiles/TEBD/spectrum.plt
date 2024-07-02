reset
complex(x,y) = x*{1,0}+y*{0,1}
target = system("echo $TARGET")
outfile = "./tex/".target."/spectrum.tex"
sfile = "snapshots/".target."/Spec/FUN/AB_left_sym.dat"
afile = "snapshots/".target."/Spec/FUN/AB_left_asym.dat"
tfile = "snapshots/".target."/Spec/FUN/AB_left_tot.dat"

nline = system("cat ".sfile." | wc -l")
array topEV[nline+0]
stats sfile u (topEV[$0+1] = abs(complex($2,$3))) nooutput
print(topEV)

set term tikz standalone size 8in,6in
set output outfile
# set title target

set xlabel 'TEBD steps'
set ylabel 'Transfer Matrix Eigenvalues'
set xrange [0:1000]
unset key

set yrange [0:1.1]
plot for [i=2:16] sfile u 0:(abs(complex(column(2*i),column(2*i+1))) / topEV[$0+1]) w l lc 1 lw 3, for [i=1:16] afile u 0:(abs(complex(column(2*i),column(2*i+1))) / topEV[$0+1]) w l lc 2 lw 2

# set yrange [0:1]
# plot for [i=1:31] file u 0:(abs(complex(column(2*i+1),column(2*i+2)))/abs(complex(column(1),column(2)))) w l

unset output
system("/Library/TeX/texbin/pdflatex -output-directory ./plots/".target." ".outfile.";")

system("echo latex")