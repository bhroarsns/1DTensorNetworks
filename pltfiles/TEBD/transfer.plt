reset
target = system("echo $TARGET")

do for [is=1:4] {
    tmfile = sprintf("./snapshots/".target."/sd%d_tm.dat", is)
    dim = system("tail -n 1 ".tmfile." | awk -F ',' '{print $1}'") + 0

    set term png size 1600,1600

    set output sprintf("./plots/".target."/tm%d.png", is)
    set title target
    set isotropic
    unset key
    set cbrange [0:1]
    set xrange [0.5:0.5+dim]
    set yrange [0.5:0.5+dim]
    plot tmfile u 1:2:3 sparse matrix=(dim,dim) origin=(1,1) with image
    unset output

    set output sprintf("./plots/".target."/logtm%d.png", is)
    set title target
    set isotropic
    unset key
    set cbrange [0:-10]
    set xrange [0.5:0.5+dim]
    set yrange [0.5:0.5+dim]
    plot tmfile u 1:2:(log($3)) sparse matrix=(dim,dim) origin=(1,1) with image
    unset output
}
