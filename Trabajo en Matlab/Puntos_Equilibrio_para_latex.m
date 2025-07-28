%Aqui despejamos directamente el valor de los puntos de equilibrio en
%funcion de c
syms k_1 mu_2 k_2 r_2 a b g_2 p_2 g_3 mu_3 c clear

%Estos vienen como resultado de los despejes para hallar los puntos de 
%equilibrio
y_2 = (-a*k_1+r_2+k_2-b*g_2*r_2)/(2*b*r_2);
y_1 = (-a*k_1+r_2-k_2-b*g_2*r_2)/(2*b*r_2);
y = y_2; %Para seleccionar la variante, correspondiente a uno de los puntos
x = k_1*y ;
z = (p_2*x*y)/(mu_3*(g_3+y));
k_1_sub = c/mu_2;
k_2_sub = sqrt(a^2*k_1^2+2*a*b*g_2*k_1*r_2-2*a*k_1*r_2+b^2*g_2^2*r_2^2+2*b*g_2*r_2^2+r_2^2);
k_2_sub = subs(k_2_sub, k_1, k_1_sub);

k_valores = [k_1;k_2];
k_valores_sub = [k_1_sub;k_2_sub];
x_s1 = subs(x, k_valores, k_valores_sub);
y_s1 = subs(y, k_valores, k_valores_sub);
z_s1 = subs(z, k_valores, k_valores_sub);

valores= [mu_2;r_2;a;b;g_2;p_2;g_3;mu_3];
sub_vals = [0.03;0.18;1;10^-9;10^5;5;10^3;10];
x_s2 = subs(x_s1,valores,sub_vals);
y_s2 = subs(y_s1,valores,sub_vals);
z_s2 = subs(z_s1,valores,sub_vals);

x_simp = simplify(x_s2);
y_simp = simplify(y_s2);
z_simp = simplify(z_s2);

disp("Simplificando")
disp(latex(x_simp))
disp(latex(y_simp))
disp(latex(z_simp))
