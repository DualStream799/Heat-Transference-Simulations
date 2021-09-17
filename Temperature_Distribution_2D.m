% Range intervals:
x_range = 0.5; % m
y_range = 0.5; % m
t_range = 10; % s
% Step intervals:
dx = 0.1;      % m
dy = 0.1;      % m
dt = 0.01;     % s
% Thermic coefficient:
alpha = 0.25;  %m²/s

% Node amounts:
x_nodes = x_range/dx + 1;
y_nodes = x_range/dy + 1;
t_nodes = t_range/dt + 1;
% 'x-axis', 'y-axis' and 'time' scales:
x_array = linspace(0, x_range, x_nodes);
y_array = linspace(0, y_range, y_nodes);
t_array = linspace(0, t_range, t_nodes);
% Temperature array:
T = zeros(x_nodes, y_nodes, t_nodes);
% Initial temperature conditions:
T(1,:,:) = 100;
T(1,1,:) = 0;
T(1, y_nodes, :) = 0;
% Iteration to calculate the temperature distribution:
for t = 2:t_nodes
    % Calculate the Temperature on the next node:
    
    for j = 2:x_nodes-1
        
        for  i = 2:y_nodes-1
            % Central Diferences:
            d2Tdx2 = ( T(i+1,j,t-1) - 2*T(i,j,t-1) + T(i-1,j,t-1) )/ dx^2;
            d2Tdy2 = ( T(i,j+1,t-1) - 2*T(i,j,t-1) + T(i,j-1,t-1) )/ dy^2;
            dTdt = alpha*dt*( d2Tdx2 + d2Tdy2 ) + T(i,j,t-1);
            T(i,j,t) = dTdt;
        end
    end
end

% Plot inital chart:
figure;
s = surf(x_array, y_array, T(:,:,1));
% s.EdgeColor = 'none';
s.FaceColor = 'interp';
view(2);
xlabel('Comprimento do eixo-x (m)');
ylabel('Comprimento do eixo-y (m)');
zlabel('Temperatura (ºC)');
axis ij;

% Plot final chart:
figure;
s = surf(x_array, y_array, T(:,:,t_nodes-1));
% s.EdgeColor = 'none';
s.FaceColor = 'interp';
view(2);
xlabel('Comprimento da barra (m)');
ylabel('Comprimento do eixo-y (m)');
zlabel('Temperatura (ºC)');
axis ij;
