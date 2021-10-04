% Range intervals:
x_range = 0.4; % m
y_range = 0.4; % m
t_range = 10; % s
% Step intervals:
dx = 0.1;      % m
dy = 0.1;      % m
dt = 0.001;     % s
% Thermic coefficient:
k = 0.23;  % W/mm.ºC
% Specific heat:
c = 897;   % J/kg.ºC
% Volumetric Density;:
d = 2.7*1e-9; % kg/mm³
% Thermic coefficient:
alpha = ( k /(d * c) )*1e-6; %m²/s
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
T(1,:,:) = 150;
T(:,y_nodes,:) = 50;

% Iteration to calculate the temperature distribution:
for t = 2:t_nodes
     
    for j = 2:y_nodes-1
        
        for i = 1:x_nodes-1
            % Calculate the Temperature on the isolated corner's nodes:
            if (i == 1)
                dTdx = 0;
                % Calculate the Central Diference on the X-axis:
                d2Tdx2 = ( T(j,i+1,t-1) - 2*T(j,i,t-1) + T(j,i,t-1) - 2*dt*dTdx )/dx^2;
                % Calculate the Central Diference on the Y-axis:
                d2Tdy2 = ( T(j+1,i,t-1) - 2*T(j,i,t-1) + T(j-1,i,t-1) )/dy^2;
                % Calculate the Temperature variation on the node:
                dTdt = alpha*dt*( d2Tdx2 + d2Tdy2);
                % Calculate the current Temperature on the node:
                T(j,i,t) = dTdt + T(j,i,t-1);
            % Calculate the Temperature on the internal nodes:
            else
                % Calculate the Central Diference on the X-axis:
                d2Tdx2 = ( T(i+1,j,t-1) - 2*T(i,j,t-1) + T(i-1,j,t-1) )/ dx^2;
                % Calculate the Central Diference on the Y-axis:
                d2Tdy2 = ( T(i,j+1,t-1) - 2*T(i,j,t-1) + T(i,j-1,t-1) )/ dy^2;
                % Calculate the Temperature variation on the node:
                dTdt = alpha*dt*( d2Tdx2 + d2Tdy2 );
                % Calculate the current Temperature on the node:
                T(i,j,t) = dTdt + T(i,j,t-1);
            end
        end
    end
end

% Plot inital chart:
figure;
s = surf(x_array, y_array, T(:,:,1));
s.EdgeColor = 'none';
s.FaceColor = 'interp';
view(2);
xlabel('Comprimento do eixo-x (m)');
ylabel('Comprimento do eixo-y (m)');
zlabel('Temperatura (ºC)');
axis ij;

% Plot final chart:
figure;
s = surf(x_array, y_array, T(:,:,t_nodes-1));
s.EdgeColor = 'none';
s.FaceColor = 'interp';
view(2);
xlabel('Comprimento da barra (m)');
ylabel('Comprimento do eixo-y (m)');
zlabel('Temperatura (ºC)');
axis ij;
