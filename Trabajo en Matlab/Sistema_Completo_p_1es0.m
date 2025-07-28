%Con este codigo estan determinados los puntos de equilibrio
syms c m_2 p_1 a r_2 m_3 p_2 g_1 g_2 g_3 b s_1 s_2 clear;
syms x y z clear;

assume(c>0)
assume(m_2>0)
assume(p_1>0)
assume(a>0)
assume(r_2>0)
assume(m_3>0)
assume(p_2>0)
assume(g_1>0)
assume(g_2>0)
assume(g_3>0)
assume(b>0)

dx = c*y - m_2*x == 0;
dy = r_2*y*(1-b*y) - (a*x*y)/(g_2+y) ==0;
dz = (p_2*x*y)/(g_3+y) - m_3*z ==0;



sols = solve([dx,dy,dz],[x,y,z]);
%disp(sols)

% 
% valores = [c,m_2,p_1,a,r_2,m_3,p_2,g_1,g_2,g_3,b];
% subs_valores = [0.5,0.03,0.1245,1,0.18,10,5,2*10^7,10^5,10^3,10^(-9)];
% dx_s = subs(dx,valores,subs_valores);
% dy_s = subs(dy,valores,subs_valores);
% dz_s = subs(dz,valores,subs_valores);
% 
% sols = solve([dx_s,dy_s,dz_s],[x,y,z]);

%sols = solve([dx,dy,dz],[x,y,z], 'Real', true);
sols_tuples = [sols.x sols.y sols.z];
disp(sols_tuples)