Graphic through gnuplot

==============================================
Start gnuplot
set GP_HOME=e:\util\gnuplot
%GP_HOME%\bin\gnuplot.exe
==============================================

==============================================
Commands to obtain the images for each example:
==============================================

plot "eq1.csv" u 1:2 title "y(t)" w lines, "eq1.csv" u 1:3 title "y'(t)" w lines

plot "eq2.csv" u 1:2 title "f(t)" w lines, "eq2.csv" u 1:3 title "theta(t)" w lines

plot "eq3.csv" u 1:2 title "u(t)" w lines, "eq3.csv" u 1:3 title "u'(t)" w lines

plot "eq4.csv" u 1:2 title "x_1(t)" w lines, "eq4.csv" u 1:3 title "x_2(t)" w lines

==============================================