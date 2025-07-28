%Este de aqui sirvio de auxilio determinando si para algun valor de c los
%puntos de equilibrio se indeterminaban
syms c

sqr = 2500000000000000*c^2 - 26997300000000*c + 72914580729;

prob = ge(sqr,0);

sols = solve(prob,c,'ReturnConditions',true, 'Real', true);
disp(sols)
disp(sols.c)
disp(sols.parameters)
disp(sols.conditions)
