reset
set term png size 640,480
target = system("echo $TARGET")
files = system("ls ./snapshots/$TARGET/Corr")
files = split(files, "\n")
nfile = |files|

set output "./plots/".target."/corr.png"
plot for [i = 1:nfile] sprintf("./snapshots/%s/Corr/%s", target, files[i]) u 1:2
unset output