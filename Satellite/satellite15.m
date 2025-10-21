%% Satellite Attitude Control with Fuzzy Logic Controller
clear all
close all
clc

%% Create results directory
resultsDir = 'Satellite_FLC_Results';
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end
timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');

%% 1. SYSTEM DEFINITION
s = tf('s');
G = 1/((s+1)*(s+9));

%% 2. LINEAR PI CONTROLLER DESIGN (for comparison)
c = 1.1;  
K = 20;  

% PI controller
Gc = K*(s + c)/s;
sys_open = Gc * G;
sys_closed = feedback(sys_open, 1);

% Step response analysis
figure(1);
step(sys_closed);
title('Linear PI Controller - Step Response');
grid on;
saveas(gcf, fullfile(resultsDir, [timestamp '_PI_Controller_StepResponse.png']));
saveas(gcf, fullfile(resultsDir, [timestamp '_PI_Controller_StepResponse.fig']));
step_info = stepinfo(sys_closed);
fprintf('Linear PI Controller:\n');
fprintf('Overshoot: %.2f%%\n', step_info.Overshoot);
fprintf('Rise Time: %.2f sec\n', step_info.RiseTime);
fprintf('Settling Time: %.2f sec\n', step_info.SettlingTime);
fprintf('Kp = %.2f, Ki = %.2f\n\n', K, K*c);

% Save PI controller results
pi_results = struct();
pi_results.Overshoot = step_info.Overshoot;
pi_results.RiseTime = step_info.RiseTime;
pi_results.SettlingTime = step_info.SettlingTime;
pi_results.Kp = K;
pi_results.Ki = K*c;

save(fullfile(resultsDir, [timestamp '_PI_Controller_Results.mat']), 'pi_results');

%% 3. OPTIMIZED FUZZY CONTROLLER DESIGN
fis = mamfis('Name', 'Satellite_FLC_Optimized');

% Input and output variables
fis = addInput(fis, [-1, 1], 'Name', 'E');
fis = addInput(fis, [-1, 1], 'Name', 'DE');
fis = addOutput(fis, [-1, 1], 'Name', 'DU');

% Membership functions for E - more aggressive response
fis = addMF(fis, 'E', 'trapmf', [-1.2, -1, -0.8, -0.5], 'Name', 'NL');
fis = addMF(fis, 'E', 'trimf', [-0.8, -0.5, -0.2], 'Name', 'NM');
fis = addMF(fis, 'E', 'trimf', [-0.5, -0.2, 0], 'Name', 'NS');
fis = addMF(fis, 'E', 'trimf', [-0.2, 0, 0.2], 'Name', 'ZR');
fis = addMF(fis, 'E', 'trimf', [0, 0.2, 0.5], 'Name', 'PS');
fis = addMF(fis, 'E', 'trimf', [0.2, 0.5, 0.8], 'Name', 'PM');
fis = addMF(fis, 'E', 'trapmf', [0.5, 0.8, 1, 1.2], 'Name', 'PL');

% Membership functions for DE
fis = addMF(fis, 'DE', 'trapmf', [-1.2, -1, -0.7, -0.4], 'Name', 'NV');
fis = addMF(fis, 'DE', 'trimf', [-0.7, -0.4, -0.2], 'Name', 'NL');
fis = addMF(fis, 'DE', 'trimf', [-0.4, -0.2, -0.1], 'Name', 'NM');
fis = addMF(fis, 'DE', 'trimf', [-0.2, -0.1, 0], 'Name', 'NS');
fis = addMF(fis, 'DE', 'trimf', [-0.1, 0, 0.1], 'Name', 'ZR');
fis = addMF(fis, 'DE', 'trimf', [0, 0.1, 0.2], 'Name', 'PS');
fis = addMF(fis, 'DE', 'trimf', [0.1, 0.2, 0.4], 'Name', 'PM');
fis = addMF(fis, 'DE', 'trimf', [0.2, 0.4, 0.7], 'Name', 'PL');
fis = addMF(fis, 'DE', 'trapmf', [0.4, 0.7, 1, 1.2], 'Name', 'PV');

% Membership functions for DU - more aggressive control
fis = addMF(fis, 'DU', 'trapmf', [-1.2, -1, -0.7, -0.4], 'Name', 'NL');
fis = addMF(fis, 'DU', 'trimf', [-0.7, -0.4, -0.2], 'Name', 'NM');
fis = addMF(fis, 'DU', 'trimf', [-0.4, -0.2, -0.1], 'Name', 'NS');
fis = addMF(fis, 'DU', 'trimf', [-0.2, -0.1, 0.1], 'Name', 'ZR');
fis = addMF(fis, 'DU', 'trimf', [-0.1, 0.1, 0.2], 'Name', 'PS');
fis = addMF(fis, 'DU', 'trimf', [0.1, 0.2, 0.4], 'Name', 'PM');
fis = addMF(fis, 'DU', 'trimf', [0.2, 0.4, 0.7], 'Name', 'PL');
fis = addMF(fis, 'DU', 'trapmf', [0.4, 0.7, 1, 1.2], 'Name', 'PV');

% Rule base 
ruleList = [...
% E   DE  DU  Weight Connection
1   1   1   1   1;   1   2   1   1   1;   1   3   1   1   1;   1   4   2   1   1;   1   5   3   1   1;
1   6   4   1   1;   1   7   5   1   1;   1   8   6   1   1;   1   9   7   1   1;

2   1   1   1   1;   2   2   1   1   1;   2   3   2   1   1;   2   4   3   1   1;   2   5   4   1   1;
2   6   5   1   1;   2   7   6   1   1;   2   8   7   1   1;   2   9   8   1   1;

3   1   1   1   1;   3   2   2   1   1;   3   3   3   1   1;   3   4   4   1   1;   3   5   5   1   1;
3   6   6   1   1;   3   7   7   1   1;   3   8   8   1   1;   3   9   8   1   1;

4   1   2   1   1;   4   2   3   1   1;   4   3   4   1   1;   4   4   5   1   1;   4   5   6   1   1;
4   6   7   1   1;   4   7   8   1   1;   4   8   8   1   1;   4   9   8   1   1;

5   1   3   1   1;   5   2   4   1   1;   5   3   5   1   1;   5   4   6   1   1;   5   5   7   1   1;
5   6   8   1   1;   5   7   8   1   1;   5   8   8   1   1;   5   9   8   1   1;

6   1   4   1   1;   6   2   5   1   1;   6   3   6   1   1;   6   4   7   1   1;   6   5   8   1   1;
6   6   8   1   1;   6   7   8   1   1;   6   8   8   1   1;   6   9   8   1   1;

7   1   5   1   1;   7   2   6   1   1;   7   3   7   1   1;   7   4   8   1   1;   7   5   8   1   1;
7   6   8   1   1;   7   7   8   1   1;   7   8   8   1   1;   7   9   8   1   1];

fis = addRule(fis, ruleList);

% Save FIS structure
writeFIS(fis, fullfile(resultsDir, [timestamp '_Satellite_FLC.fis']));

% Optimized Scaling gains
Ke = 0.15;  
Kde = 0.03;
Ki = 1.5; 

fprintf('Optimized Fuzzy Controller Gains:\n');
fprintf('Ke = %.2f, Kde = %.2f, Ki = %.2f\n\n', Ke, Kde, Ki);

%% 4. EXACT DISCRETE SIMULATION
T = 0.01;  % Sampling time
t_sim = 0:T:10;
r = 40 * ones(size(t_sim));  % Step reference

% Discretize system
sys_d = c2d(G, T, 'zoh');
[A, B, C, D] = ssdata(sys_d);

% Initialize variables
y = zeros(size(t_sim));
u = zeros(size(t_sim));
e = zeros(size(t_sim));
de = zeros(size(t_sim));
x = zeros(size(A,1), 1);  % State vector

% Simulation loop
for k = 2:length(t_sim)
    % Error calculation
    e(k) = r(k) - y(k-1);
    
    % Error derivative (filtered)
    if k > 2
        de(k) = (e(k) - e(k-1)) / T;
        de(k) = 0.9*de(k-1) + 0.1*de(k); 
    end
    
    % Normalize inputs
    e_norm = max(min(e(k) * Ke, 1), -1);
    de_norm = max(min(de(k) * Kde, 1), -1);
    
    % Fuzzy inference
    du_norm = evalfis(fis, [e_norm, de_norm]);
    
    % Control signal
    u(k) = u(k-1) + du_norm * Ki;
    
    % System update
    x = A*x + B*u(k);
    y(k) = C*x + D*u(k);
end

% Calculate performance metrics
fuzzy_overshoot = max(0, (max(y) - r(end)) / r(end) * 100);
rise_time = risetime(y, t_sim);
settling_time = settlingtime(y, t_sim, 0.02);

fprintf('Optimized Fuzzy Controller:\n');
fprintf('Overshoot: %.2f%%\n', fuzzy_overshoot);
fprintf('Rise Time: %.2f sec\n', rise_time);
fprintf('Settling Time: %.2f sec\n', settling_time);

% Rule analysis
e_test = -0.5;  % NM
de_test = 0;    % ZR
du_test = evalfis(fis, [e_test, de_test]);
fprintf('\nFor E=NM (%.2f) and DE=ZR (%.2f), output DU=%.2f\n', e_test, de_test, du_test);

% Save fuzzy controller results
fuzzy_results = struct();
fuzzy_results.Overshoot = fuzzy_overshoot;
fuzzy_results.RiseTime = rise_time;
fuzzy_results.SettlingTime = settling_time;
fuzzy_results.Ke = Ke;
fuzzy_results.Kde = Kde;
fuzzy_results.Ki = Ki;
fuzzy_results.RuleExample = struct('E', e_test, 'DE', de_test, 'DU', du_test);

save(fullfile(resultsDir, [timestamp '_Fuzzy_Controller_Results.mat']), 'fuzzy_results');

%% 5. SCENARIO 2 - RAMP RESPONSE
t2 = 0:T:16;
r2 = zeros(size(t2));

% Create ramp signal (Figure 3)
r2(t2 >= 0.4 & t2 < 5) = 60/(5-0.4) * (t2(t2 >= 0.4 & t2 < 5) - 0.4);
r2(t2 >= 5 & t2 < 8) = 60;
r2(t2 >= 8 & t2 < 16) = 60 - 60/(16-8) * (t2(t2 >= 8 & t2 < 16) - 8);

% Initialize for ramp simulation
y2 = zeros(size(t2));
u2 = zeros(size(t2));
e2 = zeros(size(t2));
de2 = zeros(size(t2));
x2 = zeros(size(A,1), 1);

% Ramp simulation
for k = 2:length(t2)
    e2(k) = r2(k) - y2(k-1);
    
    if k > 2
        de2(k) = (e2(k) - e2(k-1)) / T;
        de2(k) = 0.9*de2(k-1) + 0.1*de2(k);
    end
    
    e_norm2 = max(min(e2(k) * Ke, 1), -1);
    de_norm2 = max(min(de2(k) * Kde, 1), -1);
    
    du_norm2 = evalfis(fis, [e_norm2, de_norm2]);
    u2(k) = u2(k-1) + du_norm2 * Ki;
    
    x2 = A*x2 + B*u2(k);
    y2(k) = C*x2 + D*u2(k);
end

% Calculate tracking error (ignore initial transient)
valid_idx = t2 > 1;  % Skip first second
tracking_error = rms(r2(valid_idx) - y2(valid_idx));
fprintf('Ramp Tracking Error: %.4f degrees\n', tracking_error);

% Save ramp response results
ramp_results = struct();
ramp_results.TrackingError = tracking_error;
ramp_results.Time = t2;
ramp_results.Reference = r2;
ramp_results.Response = y2;
ramp_results.Control = u2;

save(fullfile(resultsDir, [timestamp '_Ramp_Response_Results.mat']), 'ramp_results');

%% 6. PLOTS
% Step response comparison
figure(2);
plot(t_sim, r, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Reference');
hold on;
[y_linear, t_linear] = step(sys_closed * 40, t_sim); % Scale to 40 degrees
plot(t_linear, y_linear, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Linear PI');
plot(t_sim, y, 'r-', 'LineWidth', 2, 'DisplayName', 'Fuzzy Control');
title('Step Response Comparison');
xlabel('Time (sec)');
ylabel('Angle (degrees)');
legend('Location', 'best');
grid on;
saveas(gcf, fullfile(resultsDir, [timestamp '_Step_Response_Comparison.png']));
saveas(gcf, fullfile(resultsDir, [timestamp '_Step_Response_Comparison.fig']));

% Control effort comparison
figure(3);
subplot(2,1,1);
plot(t_sim, u, 'r-', 'LineWidth', 2);
title('Fuzzy Controller - Control Effort');
xlabel('Time (sec)');
ylabel('Control Signal u(t)');
grid on;

subplot(2,1,2);
[y_step, t_step] = step(sys_closed, t_sim);
u_linear = lsim(Gc, y_step, t_step);
plot(t_step, u_linear, 'b-', 'LineWidth', 2);
title('Linear PI Controller - Control Effort');
xlabel('Time (sec)');
ylabel('Control Signal u(t)');
grid on;
saveas(gcf, fullfile(resultsDir, [timestamp '_Control_Effort_Comparison.png']));
saveas(gcf, fullfile(resultsDir, [timestamp '_Control_Effort_Comparison.fig']));

% Membership functions
figure(4);
subplot(3,1,1);
plotmf(fis, 'input', 1);
title('Input E Membership Functions');
subplot(3,1,2);
plotmf(fis, 'input', 2);
title('Input DE Membership Functions');
subplot(3,1,3);
plotmf(fis, 'output', 1);
title('Output DU Membership Functions');
saveas(gcf, fullfile(resultsDir, [timestamp '_Membership_Functions.png']));
saveas(gcf, fullfile(resultsDir, [timestamp '_Membership_Functions.fig']));

% 3D control surface
figure(5);
gensurf(fis);
title('3D Control Surface');
xlabel('Error (E)');
ylabel('Error Derivative (DE)');
zlabel('Control Change (DU)');
saveas(gcf, fullfile(resultsDir, [timestamp '_Control_Surface.png']));
saveas(gcf, fullfile(resultsDir, [timestamp '_Control_Surface.fig']));

% Ramp response
figure(6);
plot(t2, r2, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Reference');
hold on;
plot(t2, y2, 'b-', 'LineWidth', 2, 'DisplayName', 'Fuzzy Control');
title('Fuzzy Controller - Ramp Response');
xlabel('Time (sec)');
ylabel('Angle (degrees)');
legend('Location', 'best');
grid on;
saveas(gcf, fullfile(resultsDir, [timestamp '_Ramp_Response.png']));
saveas(gcf, fullfile(resultsDir, [timestamp '_Ramp_Response.fig']));

% Save all simulation data
simulation_data = struct();
simulation_data.Time = t_sim;
simulation_data.Reference = r;
simulation_data.FuzzyResponse = y;
simulation_data.FuzzyControl = u;
simulation_data.Error = e;
simulation_data.ErrorDerivative = de;

save(fullfile(resultsDir, [timestamp '_Simulation_Data.mat']), 'simulation_data');

% Save summary text file
fid = fopen(fullfile(resultsDir, [timestamp '_Summary.txt']), 'w');
fprintf(fid, 'SATELLITE ATTITUDE CONTROL SIMULATION RESULTS\n');
fprintf(fid, 'Simulation Date: %s\n\n', timestamp);
fprintf(fid, 'LINEAR PI CONTROLLER:\n');
fprintf(fid, 'Overshoot: %.2f%%\n', step_info.Overshoot);
fprintf(fid, 'Rise Time: %.2f sec\n', step_info.RiseTime);
fprintf(fid, 'Settling Time: %.2f sec\n', step_info.SettlingTime);
fprintf(fid, 'Kp = %.2f, Ki = %.2f\n\n', K, K*c);
fprintf(fid, 'FUZZY CONTROLLER:\n');
fprintf(fid, 'Overshoot: %.2f%%\n', fuzzy_overshoot);
fprintf(fid, 'Rise Time: %.2f sec\n', rise_time);
fprintf(fid, 'Settling Time: %.2f sec\n', settling_time);
fprintf(fid, 'Ke = %.2f, Kde = %.2f, Ki = %.2f\n\n', Ke, Kde, Ki);
fprintf(fid, 'RAMP RESPONSE:\n');
fprintf(fid, 'Tracking Error: %.4f degrees\n', tracking_error);
fclose(fid);

fprintf('All results saved to folder: %s\n', resultsDir);
fprintf('Simulation completed successfully!\n');

%% HELPER FUNCTIONS
function rt = risetime(y, t)
    y_ss = y(end);
    y_10 = 0.1 * y_ss;
    y_90 = 0.9 * y_ss;
    
    idx_10 = find(y >= y_10, 1);
    idx_90 = find(y >= y_90, 1);
    
    if isempty(idx_10) || isempty(idx_90)
        rt = NaN;
    else
        rt = t(idx_90) - t(idx_10);
    end
end

function st = settlingtime(y, t, tolerance)
    y_ss = y(end);
    settling_band = tolerance * y_ss;
    
    % Find when system enters and stays within settling band
    within_band = abs(y - y_ss) <= settling_band;
    
    % Find the last time it leaves the band
    st = 0;
    for i = length(t):-1:1
        if ~within_band(i)
            st = t(i);
            break;
        end
    end
end