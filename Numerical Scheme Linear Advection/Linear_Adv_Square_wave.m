clc;
clear;

% Parameters
num_points = 100;
CFL = 0.8; % Courant number
a = 1;
end_time = 10;

% Define the new initial condition
initial_condition = @(x) 0 * (x <= 0.3) + (x > 0.3 & x <= 0.7) + 0 * (x >= 0.7);
x_values = linspace(0, 1, num_points + 1);
h = 1 / num_points;
k = CFL * h / a;

% Initialize variables and plot initial conditions
time = 0;
U = zeros(2, num_points + 1);
U = [initial_condition(x_values); initial_condition(x_values)];
U_temp = U;

figure;
plot(x_values, initial_condition(x_values), 'k--'); % Initial condition plot
hold on
theoretical_plot = plot([0, a * time, a * time, 1], [initial_condition(0), initial_condition(0), initial_condition(1), initial_condition(1)], 'k-');
LW_plot = plot(x_values, U(1, :), 'bo');
LF_plot = plot(x_values, U(2, :), 'r.');
hold off
xlim([0 1]);
ylim([-0.1 1.1]);
legend('init', 'theo', 'LW', 'LF');
title_text = title(sprintf('t = %0.3f', time));
xlabel('x');
ylabel('u');

% Time evolution loop
while (time + k) < end_time
    for j = 2:num_points
        % Lax-Wendroff scheme
        U_temp(1, j) = U(1, j) - 0.5 * k / h * a * (U(1, j + 1) - U(1, j - 1)) + 0.5 * (k / h * a)^2 * (U(1, j + 1) - 2 * U(1, j) + U(1, j - 1));
        % Lax-Friedrichs scheme
        U_temp(2, j) = 0.5 * (U(2, j - 1) + U(2, j + 1)) - 0.5 * k / h * a * (U(2, j + 1) - U(2, j - 1));
    end
    
    % Apply periodic boundary conditions
    U_temp(:, 1) = U_temp(:, num_points);
    U_temp(:, num_points + 1) = U_temp(:, 2);
    
    % Update variables for the next iteration
    U = U_temp;
    time = time + k;
    set(theoretical_plot, 'XData', [0, a * time, a * time, 1]);
    set(LW_plot, 'YData', U(1, :));
    set(LF_plot, 'YData', U(2, :));
    set(title_text, 'String', sprintf('t = %0.3f', time));
    drawnow;
end
