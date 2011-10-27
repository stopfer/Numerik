set terminal png
set title "Approximation of u(t) for u'(t) = -200 * t * u(t)"
set xlabel "t"
set ylabel "u(t)"
set output "exercise_01_4.png"
plot "exercise_01_4.dat"

