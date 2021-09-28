% Range intervals:
x_range = 0.5;
t_range = 500;
% Step intervals:
dx = 0.05;
dt = 5;
% Node amounts:
x_nodes = x_range/dx + 1;
t_nodes = t_range/dt + 1;
% 'x-axis' and 'time' scales:
x_array = linspace(0,0.5,x_nodes);
t_array = linspace(0,500,t_nodes);
e_array = linspace(0,500,t_nodes);
% Thermic coefficient:
alpha = 0.0001;
% Temperature array:
T_array = zeros(t_nodes, x_nodes);
% Initial temperature conditions:
T_array(1,:) = [0, 20, 20, 20, 20, 20, 20, 20, 20, 20, 0];
% Error array:
E_array = zeros(1,t_nodes);
E_temp = zeros(1,t_nodes);

% Error array initial conditions:
E_array(1) = 0;
% Iteration to calculate the temperature distribution:
for i = 2:t_nodes
    % Set Temperature on the edges:
    T_array(i,1) = 0;
    T_array(i,x_nodes) = 0;
    % Calculate the Temperature on the next node:
    for j = 2:x_nodes-1
        T_array(i,j) = T_array(i-1,j) + alpha*(dt/dx^2)*(T_array(i-1,j+1) - 2*T_array(i-1,j) + T_array(i-1,j-1));
        % Calculate the relative error on each node:
        E_temp(j) = abs( ( T_array(i,j) - T_array(i-1,j) )/( T_array(i,j) ) );
    end
    % Geting the higher value in relative error:
    E_array(i) = max(E_temp);
end
% Temperature Plot chart:
figure;
s = surf(x_array, t_array, T_array);
s.EdgeColor = 'none';
s.FaceColor = 'interp';
view(2);
xlabel('Comprimento da barra (m)');
ylabel('Tempo (s)');
zlabel('Temperatura (ºC)');
axis ij;

% Relative Error Plot chart:
figure;
plot(t_array(2:end), E_array(2:end));
xlabel('Tempo (s)');
ylabel('Relative Error');
