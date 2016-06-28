set size square
set xlabel "# of particles per side"
set ylabel "Relative Energy Error"
set logscale y
set format y "10^{%L}"
set terminal postscript enhanced color eps "Palatino-Bold,18"
set output "p3m.eps"
plot "EnergyError.dat" u 1:2 w lp lw 1 ps 2 pt 5 notitle
