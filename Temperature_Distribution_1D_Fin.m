% Range intervals:
x_range = 0.3; % m
t_range = 10; % s
% Step intervals:
dx = 0.01;      % m
dt = 0.1;     % s
% Thermic coefficient:
k = 180;  % W/m.K
% Specific heat:
c = 896;   % J/kg.K
% Volumetric density;:
d = 2.7*1e3; % kg/m³
% Thermic coefficient:
alpha = ( k /(d * c) ); %m²/s
% Convection heat transfer coefficient:
h = 100; % W/m².K
% Initial Temperature on the base:
T_b = 100 + 273; % K
% Initial Temperature on the surroundings;
T_inf = 25 + 273; % K
% Initial Temperature on the top;
T_t = T_inf; % K
% Fin's diameter:
D = 5*1e-3; % m
% Fin's length:
L = 0.3; % m
% Fin's perimeter:
P = 2*pi*(D/2); % m
% Transversal section area:
A_tr = pi*(D/2)^2; % m²
% Time interval reference to assure numeric procedure stability condition:
dt_ref = dx^2/( alpha*( (h*P*dx^2)/(k*A_tr) + 2 ) ); % s

% Verify if selected time interval is acceptable, otherwise use reference:
if dt > dt_ref
    dt = dt_ref;
end

% Node amounts:
x_nodes = x_range/dx + 1;
t_nodes = t_range/dt + 1;
% 'x-axis', 'y-axis' and 'time' scales:
x_array = linspace(0, x_range, x_nodes);
t_array = linspace(0, t_range, t_nodes);
% Temperature array:
T = zeros(x_nodes, t_nodes);
% Initial temperature conditions:
T(:,:) = T_t;
T(1, :) = T_b;

% Error tolerance value:
tol = 1*1e-10;
E_array = zeros(1,t_nodes);
E_temp = zeros(1,t_nodes);

% Error array initial conditions:
E_array(1) = 0;
% Iteration to calculate the temperature distribution:
for l = 2:t_nodes

    % Calculate the Temperature on the next node:
    for i = 2:x_nodes-1
        dTdx_cond = ( T(i+1,l-1) - 2*T(i,l-1) + T(i-1,l-1) ) / dx^2;
        dTdx_conv = - ( (h*P)/(k*A_tr) )*( T(i,l-1) - T_inf);
        dTdt = alpha*dt*( dTdx_cond + dTdx_conv );
        T(i,l) = dTdt + T(i,l-1);
        % Calculate the relative error on each node:
        E_temp(i) = abs( ( T(i,l) - T(i,l-1) )/( T(i,l) ) );
    end
    % Geting the higher value in relative error:
    E_array(l) = max(E_temp);
%     if max(E_temp) > tol
%         break
%     end
    
end
% Temperature Plot chart:
figure;
s = surf(x_array, t_array, T.'-273);
s.EdgeColor = 'none';
s.FaceColor = 'interp';
view(2);
xlabel('Comprimento da barra (m)');
ylabel('Tempo (s)');
zlabel('Temperatura (ºC)');
axis ij;

% Temperature x Position chart:
figure;
scatter(x_array(2:x_nodes), T(1:x_nodes-1,t_nodes)-273);
xlabel('Tempo (s)');
ylabel('Temperature (ºC)');

% Relative Error Plot chart:
figure;
plot(t_array(2:end), E_array(2:end));
xlabel('Tempo (s)');
ylabel('Relative Error');
