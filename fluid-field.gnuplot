set term gif animate delay 20 size 1600,400
set output "fluid-field.gif"
set xlabel "x (m)"
set ylabel "y (m)"
set xrange [0:16]
set yrange [0:4]
set size ratio 0.25
set key off

Time = 0.1
SKIP = 1
NT = 120
dt = Time/NT

# Block dimensions from C code:
# XI = 25, XF = 57 (grid indices)
# NY/4 to 3*NY/4 for y range (grid indices)
# dx = 16.0/(328-1), dy = 4.0/(82-1)
dx = 16.0/327.0
dy = 4.0/81.0
block_x1 = 25 * dx
block_x2 = 57 * dx
block_y1 = (82/4) * dy
block_y2 = (3*82/4) * dy

# Define block 
set object 1 rect from block_x1,block_y1 to block_x2,block_y2 fc rgb "black" fs solid 1.0 front

do for [k=0:NT/SKIP] {
	t = SKIP*(k+1)*dt
	plot "fluid-field.dat" index k u 1:2:3:4 every 2:2 w vectors nohead lc rgb "blue"
  set title sprintf("2D Navier-Stokes, nu = 0.3, t = %.3f s",t)
}
