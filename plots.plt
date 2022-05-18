#plots for gnuplot
#set xlabel "K" 0,1
#set ylabel "Q error" 3.5,0
#plot "T11/ErrofcKQdt.dat" tit "dt=0.00005" w linespoints ls 4, "dtmenos/T11/ErrofcKQdt.dat" tit "dt=0.0001" w linespoints ls 2
#set term pdf

####Contour plot visto de cima
set term x11
set hidden
#unset surface  #Com ou sem superf’cie
#set cntrparam levels 20
#unset clabel
#set view 0,90
###unset view   # para desfazer
splot "Sol3DUURe.txt" w l
####End Contour plot visto de cima
#"Sol-uu-Schr-Ini.txt" w l,
#set term aqua
#set xrange [-2:3]
#set xrange [*:*]
#plot  "Sols/SolUUAbs.txt" w l, "Sols/SolVV.txt" w l, "Sols/SolWW.txt" w l


#set term pdf
#set term x11


