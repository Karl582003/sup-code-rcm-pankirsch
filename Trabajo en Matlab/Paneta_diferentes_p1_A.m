clear all
% Par√°metros fijos (seg√∫n la tabla)
mu2 = 0.03;
mu3 = 10;
p2 = 5;
a = 1;
r2 = 0.18;
b = 1e-9;
c = 0.05;  % Se usa el valor m√°ximo de c
g1 = 2e7;
g2 = 1e5;
g3 = 1e3;
s1 = 0;    % Asumimos que no hay fuente externa de c√©lulas efectoras
s2 = 0;    % Asumimos que no hay fuente externa de IL-2

% Valores de p1 a probar
%p1_values = [0, 0.001145, 0.01245, 0.1245];
p1_values = 1e1*[12.45, 1.245, 0.1245,0.01245, 0.001245,0.0001245,0 ];
% Condiciones iniciales
x0 = 1;
y0 = 1;
z0 = 1;
initial_conditions = [x0; y0; z0];

% Rango de tiempo (ajustable)
tspan = [0 100];
% ConfiguraciÛn de ode45
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);  % Tolerancias m·s estrictas

% Definir un vector de tiempo con paso controlado (ej: 1000 puntos entre 0 y 100 dÌas)
tspan = linspace(0, 100, 1000);

% Configuraci√≥n de estilos para las curvas
line_styles = {'-', '-', '-', '-.','-','-','-'};
markers = {'s','x', 'd', '^','<','+','o'};
colors = ['m','b','y', 'r', 'g', 'c','k'];
line_width = 1.0;
marker_spacing = 50;  % Espaciado para los marcadores (puntos cada 'marker_spacing' pasos)

Y1=zeros(1000,3,length(p1_values));
Y_markers=zeros(20,3,length(p1_values)); 
% Crear figura
figure;
hold on;

% Resolver el sistema para cada valor de p1
for i = 1:length(p1_values)
    p1 = p1_values(i);
    
    % Definir el sistema de ODEs
    ode_system = @(t, Y) [
        c * Y(2) - mu2 * Y(1) + (p1 * Y(1) * Y(3)) / (g1 + Y(3)) + s1;  % dx/dt
        r2 * Y(2) * (1 - b * Y(2)) - (a * Y(1) * Y(2)) / (g2 + Y(2));    % dy/dt
        (p2 * Y(1) * Y(2)) / (g3 + Y(2)) - mu3 * Y(3) + s2               % dz/dt
    ];
    
    % Soluci√≥n num√©rica con ODE45
    %[t, Y] = ode45(ode_system, tspan, initial_conditions);
    % SoluciÛn numÈrica con paso controlado
    [t, Y] = ode45(ode_system, tspan, initial_conditions, options);
   Y1(:,:,i)=Y;
end

tabla1=[];tabla2=[];tabla3=[];
for i = 1:length(p1_values)
  tabla1=[tabla1, Y1(:,1,i)];  
  tabla2=[tabla2, Y1(:,2,i)];
  tabla3=[tabla3, Y1(:,3,i)];
 end
tabla1=[t,tabla1];
tabla2=[t,tabla2];
tabla3=[t,tabla3];
for i = 1:length(p1_values)
    % A√±adir marcadores solo en algunos puntos para claridad
    t_markers = t(1:marker_spacing:end);
    Y_markers(:,:,i) = Y1(1:marker_spacing:end,:,i);
end

    % Graficar con estilos distintos para cada p1
    subplot(1, 3, 1);
    hold on
    for i = 1:length(p1_values)
    plot(t_markers, Y_markers(:, 1,i), markers{i}, 'Color', colors(i), 'MarkerSize', 6, 'LineWidth', 1);
    end
    for i = 1:length(p1_values)
    plot(t, Y1(:, 1,i), 'Color', colors(i),'markerSize',14, 'LineStyle', line_styles{i}, 'LineWidth', line_width);
     end
    hold off
    xlabel('Tiempo (dÌas)', 'FontSize',14);
    ylabel('CÈlulas efectoras (x(t))','FontSize',14);
    legend_labels = arrayfun(@(p) sprintf('p_1 = %.4f', p), p1_values, 'UniformOutput', false);
    legend(legend_labels,'FontSize',14, 'Location', 'best','NumColumns',7);
    box on
    grid on
    grid minor
    
    subplot(1, 3, 2);
    hold on
    for i = 1:length(p1_values)
    plot(t_markers, Y_markers(:, 2,i), markers{i}, 'Color', colors(i), 'MarkerSize', 6, 'LineWidth', 1);
    end
    for i = 1:length(p1_values)
    plot(t, Y1(:, 2,i), 'Color', colors(i),'markerSize',14, 'LineStyle', line_styles{i}, 'LineWidth', line_width);
    end
    hold off
   xlabel('Tiempo (dÌas)', 'FontSize',14)
    ylabel('CÈlulas tumorales (y(t))','FontSize',14);
    box on
    grid on
    grid minor
    
    subplot(1,3, 3);
    hold on
    for i = 1:length(p1_values)
    plot(t_markers, Y_markers(:, 3,i), markers{i}, 'Color', colors(i), 'MarkerSize', 6, 'LineWidth', 1);
    end
    for i = 1:length(p1_values)
    plot(t, Y1(:, 3,i), 'Color', colors(i),'markerSize',14, 'LineStyle', line_styles{i}, 'LineWidth', line_width);
    end
    hold off
    xlabel('Tiempo (dÌas)', 'FontSize',14)
    ylabel('ConcentraciÛn de IL-2 (z(t))','FontSize',14);
    box on
    grid on
    grid minor
    
    