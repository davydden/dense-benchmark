# plot parsed data from Cholesky benchmarks:
set terminal postscript eps enhanced color dashed "Helvetica" 22

############################
# line styles              #
############################
# red:
set style line 1  linetype 1 linecolor rgb "#e41a1c"  linewidth 3.000 pointtype 4 pointsize 2.0
# blue:
set style line 2  linetype 1 linecolor rgb "#377eb8"  linewidth 3.000 pointtype 6 pointsize 2.0
# green:
set style line 3  linetype 1 linecolor rgb "#4daf4a"  linewidth 3.000 pointtype 8 pointsize 2.0
# red:
set style line 11  linetype 1 dt 4 linecolor rgb "#e41a1c"  linewidth 3.000 pointtype 4 pointsize 2.0
# blue:
set style line 22  linetype 1 dt 4 linecolor rgb "#377eb8"  linewidth 3.000 pointtype 6 pointsize 2.0
# green:
set style line 33  linetype 1 dt 4 linecolor rgb "#4daf4a"  linewidth 3.000 pointtype 8 pointsize 2.0
# purple:
set style line 4  linetype 1 linecolor rgb "#984ea3"   linewidth 3.000 pointtype 9 pointsize 2.0
# orange:
set style line 5  linetype 5 dt 2 linecolor rgb "#ff7f00"   linewidth 3.000 pointtype 11 pointsize 2.0
# yellow:
set style line 6  linetype 5 dt 3 linecolor rgb "#ffff33"   linewidth 3.000 pointtype 5 pointsize 2.0
# brown
set style line 7  linetype 8 dt 4 linecolor rgb "#a65628"   linewidth 3.000 pointtype 8 pointsize 2.0
# pink
set style line 8  linetype 8 dt 5 linecolor rgb "#f781bf"   linewidth 3.000 pointtype 8 pointsize 2.0
# grey:
set style line 9  linetype 4 linecolor rgb "#999999"    linewidth 4.000 pointtype 1 pointsize 0.0
# black:
set style line 10  linetype 1 linecolor rgb "black"    linewidth 2.000 pointtype 1 pointsize 0.0

# only plot block data
extr_block(b,fname)="< awk '$2==".b."' ".fname

# a function to get a value of certain (size,block) in the given row of a given file
sb_timing(col,fname,s,b)=system("awk '$1==".s." && $2==".b." { print $" . col . "}' ".fname);

plot_lin(x,T0) = 10**(-log10(x)+log10(T0))

#--------------------------#
#                          #
# ACTUAL PLOT STARTS HERE  #
#                          #
#--------------------------#

p5_16 = "./processed_5000_32.txt"
p5_32 = "./processed_5000_16.txt"
p5_64 = "./processed_5000_64.txt"
p20_16 = "./processed_20000_16.txt"
p20_32 = "./processed_20000_32.txt"
p20_64 = "./processed_20000_64.txt"

c1 = "./cholesky_1"
c4 = "./cholesky_4"
c100 = "./cholesky_100"

lapack5 = 0 + sb_timing(5,c1,5000,16) # "0" make variable numeric !
title5 = sprintf("5000x5000 (Lapack=%.4f s)",lapack5);
title20 = sprintf("20000x20000");

sca5 =  0 + sb_timing(6,c4,5000,16)
sca20 = 0 + sb_timing(6,c4,20000,16)

# ideal scaling: T(p) = T0 * p0 / p
#
# log10 (T) = log10(T0*p0) - log10(p)
# log10 (T) = log10(T0*p0) - log2(p)/log2(10)

set autoscale xy
set xlabel "number of MPI cores"
set ylabel "wallclock [s]"
set key top right

set logscale y
set logscale x 2

set output 'cholesky.eps'
set title "strong scaling"
set size 3.0,2.0
set multiplot
set title title5
set yrange [0.3:15]
set origin 0.0,0.0
set size   1.5,2.0
plot p5_32 using (column(1)):(column(5)) axis x1y1 with lp ls 1 title 'ScaLapack (32)', \
     p5_64 using (column(1)):(column(5)) axis x1y1 with lp ls 2 title 'ScaLapack (64)', \
     p5_16 using (column(1)):(column(5)) axis x1y1 with lp ls 3 title 'ScaLapack (16)', \
     plot_lin(x,sca5*4) with l ls 10 dt 4 title 'linear scaling'
#     p5_32 using (column(1)):(column(6)) axis x1y1 with lp ls 11 title 'Elemental (32)', \
#     p5_64 using (column(1)):(column(6)) axis x1y1 with lp ls 22 title 'Elemental (64)', \

set origin 1.5,0.0
set size   1.5,2.0
set title title20
set yrange [10:250]
plot p20_32 using (column(1)):(column(5)) axis x1y1 with lp ls 1 title 'ScaLapack (32)', \
     p20_64 using (column(1)):(column(5)) axis x1y1 with lp ls 2 title 'ScaLapack (64)', \
     p20_16 using (column(1)):(column(5)) axis x1y1 with lp ls 3 title 'ScaLapack (16)', \
     plot_lin(x,sca20*4) with l ls 10 dt 4 title 'linear scaling'

unset multiplot
unset origin


#set title "Scaling with size"
#set xlabel "N^3"
#set ylabel "wallclock"
#set format x "10^{%T}"
#set output 'cholesky.eps'
#plot extr_block(32,c100) using (column(1)**3):(column(6)) axis x1y1 with lp ls 2 title 'ScaLapack (32b, 100 cores)', \
#     extr_block(32,c100) using (column(1)**3):(column(7)) axis x1y1 with lp ls 3 title 'Elemental (32b, 100 cores)', \
#     extr_block(16,c100) using (column(1)**3):(column(6)) axis x1y1 with lp ls 4 title 'ScaLapack (16b, 100 cores)', \
