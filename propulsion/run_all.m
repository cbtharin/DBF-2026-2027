%% DBF Propulsion - Run All
% Entry point for single-motor/propeller matching

clear;
clc;

%% 1. Load motor data
run('motor_db.m');

%% 2. Load propeller data
prop_file = fullfile("reference_data", "PER3_17x10E.dat");

[prop_inches, pitch_inches] = extract_dimensions(prop_file);

% Convert propeller diameter from inches to meters
prop_d = prop_inches * 0.0254;

[rpm_list, performance_data] = extract_rpm_performance(prop_file);

%% 3. Set operating conditions
voltage = 22.2;         % V
airspeed = 20;          % m/s
rho_reference = 1.212;  % kg/m^3

%% 4. Run motor-propeller solver
results = solve_motor_prop_match( ...
    motor, voltage, performance_data, ...
    rpm_list, airspeed, prop_d, rho_reference);

%% 5. Display results
fprintf('\n===== MOTOR-PROPELLER MATCH RESULTS =====\n');
disp(results);
%% Dependancy Chain
% run_all.m
% |
% |-- motor_db.m
% |
% |-- extract_dimensions.m
% |
% |-- extract_rpm_performance.m
% |
% |-- solve_motor_prop_match.m
% |
% |-- match_motor_propeller.m
% |       |
% |       |-- calculate_motor_performance.m
% |       |
% |       |-- calculate_propeller_performance.m
% |
% |-- calculate_motor_performance.m
% |
% |-- calculate_propeller_performance.m