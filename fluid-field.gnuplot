set term gif animate delay 8 size 900,350
set output "fluid-field.gif"
set xlabel "x (m)"
set ylabel "y (m)"
set xrange [0:6]
set yrange [0:2]
set size ratio 0.33
set key off

# Obstacle - vertical rectangle
# Center (1.26, 1.0), width 0.15, height 0.7
set object 1 rect from 1.185,0.65 to 1.335,1.35 fc rgb "black" fs solid

Time = 4.0
SKIP = 20
NT = 4000
dt = Time/NT

do for [k=0:NT/SKIP] {
    t = SKIP*k*dt
    set title sprintf("Flow Past Obstacle, nu = 0.08, t = %.2f s", t)
    plot "fluid-field.dat" index k u 1:2:3:4 every 2:2 w vectors head filled lc "blue"
}
