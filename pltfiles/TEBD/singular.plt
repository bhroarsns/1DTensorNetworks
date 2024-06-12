target = system("echo $TARGET")
sitetype = 4
# titlestr = system("echo ".target." | sed 's|\/|, |g' | sed 's|1_2|1\/2|' | sed 's|_|\ |g'")
filename(is, it) = sprintf("./snapshots/%s/%c%d.dat", target, 65+it, is)

count = 0
do for [is=1:sitetype] {
    do for [it=0:1] {
        valarr = split(system(sprintf("tail -n 1 %s", filename(is, it))), ", ")
        count = count + |valarr|
    }
}

array NA[count]
array SA[count]
array RA[count]
array TA[count]

iel = 1
do for [is=1:sitetype] {
    do for [it=0:1] {
        valarr = split(system(sprintf("tail -n 1 %s", filename(is, it))), ", ")
        do for [id=1:|valarr|] {
            TA[iel] = |valarr|
            NA[iel] = id
            SA[iel] = is * int((it - 0.5) * 2)
            RA[iel] = valarr[id] + 0
            iel = iel + 1
        }
    }
}

set term tikz standalone size 6in,6in
# set term png size 1600,1600
set title target noenhanced

unset key
set grid
if (sitetype == 2) {
    set xrange [-2.5:2.5]
    set xtics ("↓" -2, "↑" -1, "↑" 1, "↓" 2)
} else if (sitetype == 3) {
    set xrange [-3.5:3.5]
    set xtics ("-" -3, "0" -2, "+" -1, "+" 1, "0" 2, "-" 3)
} else if (sitetype == 4) {
    set xrange [-4.5:4.5]
    set xtics ("UpDn" -4, "Dn" -3, "Up" -2, "Emp" -1, "Emp" 1, "Up" 2, "Dn" 3, "UpDn" 4)
}

outfile = "./tex/".target."/singular.tex"
set output outfile
do for [i=1:count] {
    if (NA[i] % 2 == 1) {
        set label sprintf("%d", NA[i]) at SA[i]-0.15-0.35*((NA[i]+0.0)/TA[i]),RA[i]+(0.5-((NA[i]-1)/2)%2)*0.01 left font ",8"
    } else {
        set label sprintf("%d", NA[i]) at SA[i]+0.15+0.35*((NA[i]+0.0)/TA[i]),RA[i]+(0.5-((NA[i]-1)/2)%2)*0.01 right font ",8"
    }
}
plot NA u (SA[$1]):(RA[$1]) ps 2 pt 6 lc 1, NA u (SA[$1]):(RA[$1]) ps 2 pt 1 lc 1

unset output
system("/Library/TeX/texbin/pdflatex -output-directory ./plots/".target." ".outfile.";")
