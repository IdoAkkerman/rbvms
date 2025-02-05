set multiplot layout 1,2
CD = 1.411
p2p = 0.0203
CL = 0.7267
plot[0:][1.38:1.44] "output_000000.dat" u 2:(-2*$9) t "C_D" w l, CD - p2p/2.0, CD, CD + p2p/2.0
plot[0:][-0.5:0.5] "output_000000.dat" u 2:(2*$10) t "C_L" w l, -CL/2.0, CL/2.0
pause -1

