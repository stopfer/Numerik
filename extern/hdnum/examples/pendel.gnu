set terminal postscript eps color 18
set output "pendel.eps"
set title "Pendel bei kleiner Auslenkung"
plot "pendel.dat" title "" w l lw 3
