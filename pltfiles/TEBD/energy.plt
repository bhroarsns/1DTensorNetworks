model = system("echo $MODEL")
target = model."/iTEBD/mpslen=2/D=16"
datfiles = system("ls results/".target."/*/*/energy.dat results/".target."/*/*/*/energy.dat")
outfile = "tex/".target."/energy.tex"

titlefn(dat) = sprintf("%s %s", system("echo ".dat." | grep -oE \"seed=[0-9]+\" | sed \"s/seed=//g\""), system("echo ".dat." | grep -oE 'initΔτ\=0\.1.*\/energy.dat' | sed \"s/initΔτ=0.1//g\" | sed \"s/\\\/energy.dat//g\" | sed \"s/\\\///g\""))

set term tikz standalone size 8in,6in
set output outfile
set title target noenhanced

set xlabel 'TEBD steps'
set ylabel 'Energy Density $E_{GS}$'
# set xrange [0:500]
set yrange [floor(exac*10) / 10.0:exac+0.1]

plot for [datfile in datfiles] datfile u 1:3 w l title titlefn(datfile), exac w l title "Exact"

set xrange [0:1]
unset yrange
set logscale y
set xlabel 'inverse of accumulated $\tau$'
set ylabel 'Relative Error Energy Density $\Delta E_{GS}$'
# set xrange [0:500]
# set yrange [exac-0.1:exac+0.1]
# set logscale y

plot for [datfile in datfiles] datfile u (1.0/$2):(abs($3-exac)/abs(exac)) w l title titlefn(datfile)


# set xlabel 'TEBD steps'
# set ylabel 'Relative Error of Energy Density $\Delta E_{GS}$'
# set logscale y
# set xrange [0:500]

# plot datfile u 1:(abs($2-exac)/abs(exac)) w l title "TEBD result"

# set xlabel 'Temperature $T$'
# set ylabel 'Relative Error of Energy Density $\Delta E_{GS}$'
# unset logscale
# set xrange [0:1.0]
# set yrange [0:]
# plot datfile u (1/($1*0.1)):(abs($2-exac)/abs(exac)) w l title "TEBD result"

# set xrange [0:0.1]
# plot datfile u (1/($1*0.1)):(abs($2-exac)/abs(exac)) w l title "TEBD result"

unset output
system("/Library/TeX/texbin/pdflatex -output-directory ./plots/".target." ".outfile.";")

system("echo latex")