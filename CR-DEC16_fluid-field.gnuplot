set term gif animate delay 10 size 1600,1200
set output "fluid-field.gif"
set xlabel "x (m)"
set ylabel "y (m)"
set xrange [0:30]
set yrange [0:20]
set size ratio 0.3
set key off

dx = 30.0/299.0
dy = 20.0/199.0

block_x1 = 20*dx
block_x2 = 25*dx
block_y1 = 70*dy
block_y2 = 130*dy
 
Time = 0.1
SKIP = 4
NT = 500
dt = Time/NT

set object 1 rect from block_x1,block_y1 to block_x2,block_y2 fc rgb "gray" fs solid 1.0 front

do for [k=0:NT/SKIP] {
	t = SKIP*(k+1)*dt
	plot "fluid-field.dat" index k u 1:2:3:4 every 2:2 w vectors nohead lc rgb"blue" 
  set title sprintf("2D Navier-Stokes, nu = 0.3, t = %.3f s",t)
}
