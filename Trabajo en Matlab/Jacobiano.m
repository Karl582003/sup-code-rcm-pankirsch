%En este script calculamos las evaluaciones del Jacobiano en los puntos
syms x y z clear
syms c mu_2 p_1 a r_2 mu_3 p_2 g_1 g_2 g_3 b s_1 s_2 clear

%Este resultado fue obtenido del script de los puntos de equilibrio
%Para luego pasarlo a un .mat para cargarlo en otros scripts
load('Puntos_Equilibrio_Coords.mat')

dx = c*y - mu_2*x;
dy = r_2*y*(1-b*y) - (a*x*y)/(g_2+y);
dz = (p_2*x*y)/(g_3+y) - mu_3*z;

valores = [c,mu_2,p_1,a,r_2,mu_3,p_2,g_1,g_2,g_3,b];
subs_valores = [0.025,0.03,0.1245,1,0.18,10,5,2*10^7,10^5,10^3,10^(-9)];

x_1 = subs(x_eq1,valores,subs_valores);
x_2 = subs(x_eq2,valores,subs_valores);
y_1 = subs(y_eq1,valores,subs_valores);
y_2 = subs(y_eq2,valores,subs_valores);
z_1 = subs(z_eq1,valores,subs_valores);
z_2 = subs(z_eq2,valores,subs_valores);


jacob = jacobian([dx,dy,dz],[x,y,z]);
% disp('Jacobiano con parametros')
% disp(latex(jacob))


jacob_subst = subs(jacob,valores,subs_valores);
% disp('Jacobiano sustituido, tomando c=0.025')
% disp(latex(jacob_subst))

jacob_1 = subs(jacob_subst, [x,y,z], [x_1,y_1,z_1]);
jacob_1 = double(jacob_1);
jacob_2 = subs(jacob_subst, [x,y,z], [x_2,y_2,z_2]);
jacob_2 = double(jacob_2);
disp('Jacobiano 1')
disp(jacob_1)
disp('Jacobiano 2')
disp(jacob_2)

P_1=double([x_1,y_1,z_1]);
P_2=double([x_2,y_2,z_2]);
disp("P_1")
disp(P_1)
disp("P_2")
disp(P_2)


