set terminal postscript eps color 18
set output "pendelnumerisch.eps"
set title "Numerische Loesung des Pendels"
plot "pendelnumerisch.dat" title "" w l lw 3
