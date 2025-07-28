%Este de aqui sirvio de auxilio determinando si para algun valor de c los
%puntos de equilibrio se indeterminaban
syms c

sqr = 10000/9 * c^2 - (29997)/(2500)*c+(8101620081)/(2.5*10^11);
den = (81*(2500000000*c + 50*(2500000000000000*c^2 - 26997300000000*c + 72914580729)^(1/2) - 13498677));

prob = den==0;
sols = solve(prob,c);
disp(sols)

c_t = -10:0.01:10;
den = subs(den,c,c_t);
plot(c,den)