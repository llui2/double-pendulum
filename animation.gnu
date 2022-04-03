set term gif size 1024,1024 animate  delay 2
set output "dpendulum.gif"

set encoding utf8

set tmargin at screen 0.99
set bmargin at screen 0.01
set rmargin at screen 0.99
set lmargin at screen 0.01

set key  top right font ",15"
set key spacing 1.5

set title "" font "Helvetica Bold,15"

set xtics font ", 15"
set ytics font ", 15"

set xrange [-2.1:2.1]
set yrange [-2.2:2.2]


set style line 1 lc rgb 'dark-green' pt 7 ps 3 lt 1 lw 2



set format x ""
set format y ""

set xlabel "" font ",15" offset 0
set ylabel "" font ",15" offset 0


do for [i=1:1000] {   plot  "dades.dat" every::i-800::i u 4:5 t '' with line lt 1 lw 1 lc rgb 'green',\
                            "dades.dat" every::i::i u 2:3 t '' with points ls 1,\
                            "dades.dat" every::i::i u 4:5 t '' with points ls 1,\
                            "dades.dat" every::i::i u 2:3:(0-$2):(0-$3) t '' with vectors nohead lw 2 lc rgb 'dark-green',\
                            "dades.dat" every::i::i u 2:3:($4-$2):($5-$3) t '' with vectors nohead lw 2 lc rgb 'dark-green',\
                            for [j=1:8] "dades.dat" every::i-j::i u 4:5 t '' with line lt 1 lw 12-j lc rgb 'dark-green'}