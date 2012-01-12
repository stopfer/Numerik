set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
#set xr [-3:3]
#set yr [-0.1:1.1]

set terminal pdf
set title "Approximation of reindeer population u(t)"
set xlabel "t"
set ylabel "u(t)"
set output "reindeerpopulation.pdf"
plot "reindeerpopulation.dat" with lines

