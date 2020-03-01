set terminal pngcairo size 1000, 1000

set output 'mesh.png'

plot 'meshLines.csv' using 1:2 notitle w l lw 2 lc 'black'
