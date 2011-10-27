set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set xr [-3:3]
set yr [0:1]

set terminal postscript enhanced color
set title "Approximation of u(t) for u'(t) = -200 * t * u(t)"
set xlabel "t"
set ylabel "u(t)"
set output "exercise_01_4.ps"
plot "exercise_01_4.dat" with lines

