set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
#set xr [-3:3]
#set yr [-0.1:1.1]

set terminal postscript enhanced color
set title "Approximation of u(t) for u'(t) = -t / u(t)"
set xlabel "t"
set ylabel "u(t)"
set output "exercise_03_3.ps"
plot "exercise_03_3.dat" with lines

