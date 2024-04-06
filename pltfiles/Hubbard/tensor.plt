set term gif animate delay 1 size 640,480
set output "./plots/F_XXX/S=1_2/iTEBD/mpslen=2/D=16/seed=10/tensor.gif"
sta(i) = sprintf("snapshots/F_XXX/S=1_2/iTEBD/mpslen=2/D=16/seed=10/A/%d.dat", i)
stb(i) = sprintf("snapshots/F_XXX/S=1_2/iTEBD/mpslen=2/D=16/seed=10/B/%d.dat", i)
bwab(i) = sprintf("snapshots/F_XXX/S=1_2/iTEBD/mpslen=2/D=16/seed=10/AB/%d.dat", i)
bwba(i) = sprintf("snapshots/F_XXX/S=1_2/iTEBD/mpslen=2/D=16/seed=10/BA/%d.dat", i)

do for [i=0:500] {
    set isotropic
    unset key
    chi=16

    set xrange [0.5:(2*chi+3.5)]
    set yrange [0.5:(2*chi+3.5)]
    stsiz=chi*chi*2
    bwsiz=chi
    ccord=2+chi

    array XC[stsiz*2+bwsiz*4]
    array YC[stsiz*2+bwsiz*4]
    array ZC[stsiz*2+bwsiz*4]
    stats sta(i)  using (XC[1+$0]                 = ccord-(1+$2), YC[1+$0]                 = ccord-(1+$3)*(3-2*$1), ZC[1+$0]                 = $6) nooutput
    stats stb(i)  using (XC[1+$0+1*stsiz]         = ccord+(1+$3), YC[1+$0+1*stsiz]         = ccord-(1+$2)*(3-2*$1), ZC[1+$0+1*stsiz]         = $6) nooutput
    stats bwab(i) using (XC[1+$0+2*stsiz]         = ccord,        YC[1+$0+2*stsiz]         = ccord-(1+$1),          ZC[1+$0+2*stsiz]         = 10*$2) nooutput
    stats bwab(i) using (XC[1+$0+2*stsiz+1*bwsiz] = ccord,        YC[1+$0+2*stsiz+1*bwsiz] = ccord+(1+$1),          ZC[1+$0+2*stsiz+1*bwsiz] = 10*$2) nooutput
    stats bwba(i) using (XC[1+$0+2*stsiz+2*bwsiz] = ccord-(1+$1), YC[1+$0+2*stsiz+2*bwsiz] = ccord,                 ZC[1+$0+2*stsiz+2*bwsiz] = 10*$2) nooutput
    stats bwba(i) using (XC[1+$0+2*stsiz+3*bwsiz] = ccord+(1+$1), YC[1+$0+2*stsiz+3*bwsiz] = ccord,                 ZC[1+$0+2*stsiz+3*bwsiz] = 10*$2) nooutput

    set cbrange [0:20]
    set title sprintf("TEBD step: %d", i)
    plot ZC u (XC[$1]):(YC[$1]):(ZC[$1]) sparse matrix=(2*chi+3,2*chi+3) origin=(1,1) with image
    # pause 0.01
}
# pause -1

system("echo gif")