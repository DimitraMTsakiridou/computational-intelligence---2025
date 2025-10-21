%% Car Control C Simulation with Fuzzy Logic
close all;
clear;
clc;

%% Variables declaration - according to problem specifications
% Initial position
x_init = 3.8;      
y_init = 0.5;      
% Car velocity
velocity = 0.05;
% Initial velocity angles (in degrees) 
theta_deg = [0, 45, -45];
% Obstacle Bounds
obstacleBounds = [[5 5 6 6 7 7 10 10];
                  [0 1 1 2 2 3 3 0]];
% Destination coordinates
x_final = 10;
y_final = 3.2;
% Final point error tolerance
errorTolerance = 0.03;
% Set bounds for runtime (prevent infinite loops)
x_bounds = [-1 (x_final + 0.3)];
y_bounds = [-1 (y_final + 0.3)];

% Create results directory if it doesn't exist
if ~exist('results', 'dir')
    mkdir('results');
end

%% (Controller with initial membership functions)
fprintf('=== Simulation with Initial Membership Functions ===\n');
carFLC = readfis("carFLC.fis");
showrule(carFLC)

% Initialize results structure
results_initial = struct();

%% ==============================================
for i = 1:length(theta_deg)
    fprintf("For initial velocity angle = %d [deg]\n", theta_deg(i))

    % Initialize variables
    x = x_init;
    y = y_init;
    theta = theta_deg(i);
    path_length = 0;
    iteration = 0;
    
    % Store initial conditions
    results_initial(i).theta_initial = theta_deg(i);
    results_initial(i).x_path = x;
    results_initial(i).y_path = y;
    results_initial(i).dv_history = [];
    results_initial(i).dh_history = [];
    results_initial(i).theta_history = theta;
    results_initial(i).dtheta_history = [];

    finalPointError = pdist([[x_final y_final];[x_init y_init]]);
    
    while(finalPointError > errorTolerance && iteration < 1000) % If near final destination, stop
        iteration = iteration + 1;
        
        % Get sensor distances from obstacle
        [dh, dv] = getSensorDistances(x(end), y(end));
        if (isnan(dh) || isnan(dv))
            fprintf("Car hit the obstacle!\n\n");
            results_initial(i).hit_obstacle = true;
            break;
        end

        % Use FLC to calculate dtheta
        dtheta = evalfis(carFLC, [dv, dh, theta]);
        
        % Calculate new velocity's angle
        theta = theta + dtheta;
        
        % Find new points using velocity's new angle
        x_new = x(end) + cosd(theta)*velocity;
        y_new = y(end) + sind(theta)*velocity;
        
        % Update path length
        path_length = path_length + sqrt((x_new - x(end))^2 + (y_new - y(end))^2);

        % Update car's path
        x = [x x_new];
        y = [y y_new];
        
        % Store data for analysis
        results_initial(i).dv_history = [results_initial(i).dv_history, dv];
        results_initial(i).dh_history = [results_initial(i).dh_history, dh];
        results_initial(i).theta_history = [results_initial(i).theta_history, theta];
        results_initial(i).dtheta_history = [results_initial(i).dtheta_history, dtheta];

        % Check if car hit the obstacle or out of bounds
        if (x(end) < x_bounds(1) || x(end) > x_bounds(2) || y(end) < y_bounds(1) || y(end) > y_bounds(2))
            fprintf("OUT OF BOUNDS. Runtime stopped...\n");
            results_initial(i).out_of_bounds = true;
            break;
        end

        % Update final point error for check
        finalPointError = pdist([[x_final y_final];[x(end) y(end)]]);
    end
    
    % Store results
    results_initial(i).x_path = x;
    results_initial(i).y_path = y;
    results_initial(i).final_error = finalPointError;
    results_initial(i).path_length = path_length;
    results_initial(i).iterations = iteration;
    results_initial(i).success = (finalPointError <= errorTolerance);
    
    fprintf("Final point error = %f [m]\n", finalPointError)
    fprintf("Path length = %f [m]\n", path_length)
    fprintf("Iterations = %d\n", iteration)
    fprintf("Success: %d\n", results_initial(i).success)
    fprintf("===========\n")

    % Plot the simulation environment
    figure;
    hold on;
    plot(polyshape(obstacleBounds(1,:), obstacleBounds(2,:)), 'FaceColor', '#808080')   % Bounds plot
    scatter(x_init, y_init, 100, 'o', 'filled', 'MarkerFaceColor','b', 'MarkerEdgeColor', 'k')  % Initial point
    scatter(x_final, y_final, 100, 's', 'filled', 'MarkerFaceColor','r', 'MarkerEdgeColor', 'k') % Destination
    plot(x, y, 'b-', "LineWidth", 1.5);                                 % Car path
    xlim([0 10.5])
    ylim([0 4])
    grid on;
    xlabel('x position [m]');
    ylabel('y position [m]');
    title("Car Control Simulation (Initial MFs)", 'FontSize', 14)
    subtitle(sprintf("Initial angle: %d°", theta_deg(i)), 'FontSize', 12)
    legend("Obstacle", "Initial Position", "Desired Position", "Car's Path", 'Location', 'northwest')
    
    % Save figure
    saveas(gcf, sprintf('results/initial_path_theta_%d.png', theta_deg(i)));
    savefig(sprintf('results/initial_path_theta_%d.fig', theta_deg(i)));
end

% Save initial results
save('results/results_initial.mat', 'results_initial');

%% (Controller with modified membership functions)
fprintf('\n\n=== Simulation with Modified Membership Functions ===\n');
carFLC_modified = readfis("carFLC_modified.fis");
showrule(carFLC_modified)

% Initialize results structure
results_modified = struct();

%% ==============================================
for i = 1:length(theta_deg)
    fprintf("For initial velocity angle = %d [deg]\n", theta_deg(i))

    % Initialize variables
    x = x_init;
    y = y_init;
    theta = theta_deg(i);
    path_length = 0;
    iteration = 0;
    
    % Store initial conditions
    results_modified(i).theta_initial = theta_deg(i);
    results_modified(i).x_path = x;
    results_modified(i).y_path = y;
    results_modified(i).dv_history = [];
    results_modified(i).dh_history = [];
    results_modified(i).theta_history = theta;
    results_modified(i).dtheta_history = [];

    finalPointError = pdist([[x_final y_final];[x_init y_init]]);
    
    while(finalPointError > errorTolerance && iteration < 1000) % If near final destination, stop
        iteration = iteration + 1;
        
        % Get sensor distances from obstacle
        [dh, dv] = getSensorDistances(x(end), y(end));
        if(isnan(dh) || isnan(dv))
            fprintf("Car hit the obstacle!\n\n");
            results_modified(i).hit_obstacle = true;
            break;
        end

        % Use FLC to calculate dtheta
        dtheta = evalfis(carFLC_modified, [dv, dh, theta]);
        
        % Calculate new velocity's angle
        theta = theta + dtheta;
        
        % Find new points using velocity's new angle
        x_new = x(end) + cosd(theta)*velocity;
        y_new = y(end) + sind(theta)*velocity;
        
        % Update path length
        path_length = path_length + sqrt((x_new - x(end))^2 + (y_new - y(end))^2);

        % Update car's path
        x = [x x_new];
        y = [y y_new];
        
        % Store data for analysis
        results_modified(i).dv_history = [results_modified(i).dv_history, dv];
        results_modified(i).dh_history = [results_modified(i).dh_history, dh];
        results_modified(i).theta_history = [results_modified(i).theta_history, theta];
        results_modified(i).dtheta_history = [results_modified(i).dtheta_history, dtheta];

        % Check if car hit the obstacle or out of bounds
        if ((x(end) < x_bounds(1) || x(end) > x_bounds(2)) || (y(end) < y_bounds(1) || y(end) > y_bounds(2)))
            fprintf("OUT OF BOUNDS. Runtime stopped...\n");
            results_modified(i).out_of_bounds = true;
            break;
        end

        % Update final point error for check
        finalPointError = pdist([[x_final y_final];[x(end) y(end)]]);
    end
    
    % Store results
    results_modified(i).x_path = x;
    results_modified(i).y_path = y;
    results_modified(i).final_error = finalPointError;
    results_modified(i).path_length = path_length;
    results_modified(i).iterations = iteration;
    results_modified(i).success = (finalPointError <= errorTolerance);
    
    fprintf("Final point error = %f [m]\n", finalPointError)
    fprintf("Path length = %f [m]\n", path_length)
    fprintf("Iterations = %d\n", iteration)
    fprintf("Success: %d\n", results_modified(i).success)
    fprintf("===========\n")

    % Plot the simulation environment
    figure;
    hold on;
    plot(polyshape(obstacleBounds(1,:), obstacleBounds(2,:)), 'FaceColor', '#808080')   % Bounds plot
    scatter(x_init, y_init, 100, 'o', 'filled', 'MarkerFaceColor','b', 'MarkerEdgeColor', 'k')  % Initial point
    scatter(x_final, y_final, 100, 's', 'filled', 'MarkerFaceColor','r', 'MarkerEdgeColor', 'k') % Destination
    plot(x, y, 'b-', "LineWidth", 1.5);                                 % Car path
    xlim([0 10.5])
    ylim([0 4])
    grid on;
    xlabel('x position [m]');
    ylabel('y position [m]');
    title("Car Control Simulation (Modified MFs)", 'FontSize', 14)
    subtitle(sprintf("Initial angle: %d°", theta_deg(i)), 'FontSize', 12)
    legend("Obstacle", "Initial Position", "Desired Position", "Car's Path", 'Location', 'northwest')
    
    % Save figure
    saveas(gcf, sprintf('results/modified_path_theta_%d.png', theta_deg(i)));
    savefig(sprintf('results/modified_path_theta_%d.fig', theta_deg(i)));
end

% Save modified results
save('results/results_modified.mat', 'results_modified');

%% Generate comparison plots
fprintf('\n\n=== Generating Comparison Plots ===\n');

% Plot all paths together for comparison
figure;
hold on;
plot(polyshape(obstacleBounds(1,:), obstacleBounds(2,:)), 'FaceColor', '#808080')   % Bounds plot
scatter(x_init, y_init, 100, 'o', 'filled', 'MarkerFaceColor','b', 'MarkerEdgeColor', 'k')  % Initial point
scatter(x_final, y_final, 100, 's', 'filled', 'MarkerFaceColor','r', 'MarkerEdgeColor', 'k') % Destination

colors = {'b', 'g', 'm'};
line_styles = {'-', '--'};
labels = {};

% Plot initial MF paths
for i = 1:length(theta_deg)
    plot(results_initial(i).x_path, results_initial(i).y_path, ...
        'Color', colors{i}, 'LineStyle', line_styles{1}, 'LineWidth', 1.5);
    labels{end+1} = sprintf('Initial MFs, θ=%d°', theta_deg(i));
end

% Plot modified MF paths
for i = 1:length(theta_deg)
    plot(results_modified(i).x_path, results_modified(i).y_path, ...
        'Color', colors{i}, 'LineStyle', line_styles{2}, 'LineWidth', 1.5);
    labels{end+1} = sprintf('Modified MFs, θ=%d°', theta_deg(i));
end

xlim([0 10.5])
ylim([0 4])
grid on;
xlabel('x position [m]');
ylabel('y position [m]');
title("Comparison of Car Paths", 'FontSize', 14)
legend([{'Obstacle', 'Start', 'Destination'}, labels], 'Location', 'northwest')

% Save comparison figure
saveas(gcf, 'results/comparison_all_paths.png');
savefig('results/comparison_all_paths.fig');

%% Display results summary
fprintf('\n=== RESULTS SUMMARY ===\n');
fprintf('Initial Membership Functions:\n');
for i = 1:length(theta_deg)
    fprintf('  θ=%d°: Error=%.4fm, Path Length=%.2fm, Success=%d\n', ...
            theta_deg(i), results_initial(i).final_error, ...
            results_initial(i).path_length, results_initial(i).success);
end

fprintf('\nModified Membership Functions:\n');
for i = 1:length(theta_deg)
    fprintf('  θ=%d°: Error=%.4fm, Path Length=%.2fm, Success=%d\n', ...
            theta_deg(i), results_modified(i).final_error, ...
            results_modified(i).path_length, results_modified(i).success);
end

fprintf('\nSimulation complete. Results saved in ''results'' folder.\n');