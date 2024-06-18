reset
set term png size 1600,1600
target = system("echo $TARGET")
files = system("ls ./snapshots/$TARGET/Corr")
files = split(files, "\n")
nfile = |files|

set output "./plots/".target."/corr.png"
# set key outside
unset key
do for [i = 1:nfile] {
    pos = split(system("tail -n 1 ./snapshots/".target."/Corr/".files[i]), ", ")
    set arrow from pos[1]+10.0,pos[2]+(i*1.0)/(nfile*1.0)*0.1 to pos[1]+0.0,pos[2]+0.0
    set label system("echo ".files[i]." | sed s/.dat//g") at pos[1]+10.0,pos[2]+(i*1.0)/(nfile*1.0)*0.1
}

plot for [i = 1:nfile] "./snapshots/".target."/Corr/".files[i] u 1:2 w l lw 3

unset output