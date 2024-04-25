target = "AF_XXX/S=1/iTEBD/mpslen=2/D=16/seed=10"

toutfile = "./plots/".target."/tm.png"
tmfile = "./snapshots/".target."/tm.dat"
bonddim = system("awk -F \", \" 'NR==1{print $1}' ".tmfile)
array TXC[bonddim * bonddim * bonddim * bonddim]
array TYC[bonddim * bonddim * bonddim * bonddim]
array TZC[bonddim * bonddim * bonddim * bonddim]
array TAC[bonddim * bonddim * bonddim * bonddim]
stats tmfile every ::1 u (TXC[1+$0] = $5, TYC[1+$0] = $6, TZC[1+$0] = $9, TAC[1+$0]=$10) nooutput

set term png size 1600,1600
set output toutfile
set title target
set isotropic
unset key
# set cbrange [0:1]

set xrange [0.5:0.5+bonddim*bonddim]
set yrange [0.5:0.5+bonddim*bonddim]
plot TZC u (TXC[$1]):(TYC[$1]):(TZC[$1]) sparse matrix=(bonddim*bonddim,bonddim*bonddim) origin=(1,1) with image

toutfile = "./plots/".target."/tmangle.png"
set term png size 1600,1600
set output toutfile
set title target
set isotropic
unset key
# set cbrange [0:1]

set xrange [0.5:0.5+bonddim*bonddim]
set yrange [0.5:0.5+bonddim*bonddim]
plot TAC u (TXC[$1]):(TYC[$1]):(TAC[$1]) sparse matrix=(bonddim*bonddim,bonddim*bonddim) origin=(1,1) with image

soutfile = "./plots/".target."/symtm.png"
sfile = "./snapshots/".target."/symtm.dat"
symdim = system("awk -F \", \" 'NR==1{print $1}' ".sfile)
array SXC[symdim * symdim]
array SYC[symdim * symdim]
array SZC[symdim * symdim]
stats sfile  every ::1 u (SXC[1+$0] = $1, SYC[1+$0] = $2, SZC[1+$0] = $5) nooutput

afile = "./snapshots/".target."/asymtm.dat"
asymdim = system("awk -F \", \" 'NR==1{print $1}' ".afile)
array AXC[asymdim * asymdim]
array AYC[asymdim * asymdim]
array AZC[asymdim * asymdim]
stats afile  every ::1 u (AXC[1+$0] = symdim+$1, AYC[1+$0] = bonddim+$2, AZC[1+$0] = $5) nooutput

set term png size 1600,1600
set output soutfile
set isotropic
unset key
# set cbrange [0:1]

set xrange [0.5:0.5+symdim+asymdim]
set yrange [0.5:0.5+symdim]
plot SZC u (SXC[$1]):(SYC[$1]):(SZC[$1]) sparse matrix=(symdim,symdim) origin=(1,1) with image, AZC u (AXC[$1]):(AYC[$1]):(AZC[$1]) sparse matrix=(asymdim,asymdim) origin=(symdim+1,bonddim+1) with image

unset output