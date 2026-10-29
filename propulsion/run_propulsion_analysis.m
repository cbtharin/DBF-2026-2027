%% DBF 2026-2027 Propulsion Analysis
% Master script for motor and propeller evaluation

clear;
clc;

%% load motor data
motor_db;

%% 1. Load propeller data

prop_file = fullfile("reference_data", "PER3_17x10E.dat");

% Extract propeller dimensions
[prop_inches, pitch_inches] = extract_dimensions(prop_file);

% Convert diameter to meters
prop_d = prop_inches * 0.0254;

% Extract APC performance data
[rpm_list, performance_data] = extract_rpm_performance(prop_file);

% Reference air density
rho_reference = 1.212;

fprintf('Loaded APC %gx%gE propeller data.\n', ...
    prop_inches, pitch_inches); %need to fix it outputs d in meters in sentence


%% 2. Define test motor

motor.name = "Scorpion SII-4025-520";
motor.Kv = 520;
motor.Rm = 0.016;
motor.I0 = 1.40;
motor.Imax = 100;
motor.Pmax = 2000;
motor.mass = 0.353;


%% 3. Test propeller performance

airspeed = 20;  % m/s

rpm_test = [7000, 7450];

for i = 1:length(rpm_test)

    rpm = rpm_test(i);

    [thrust, shaft_power, torque] = calculate_propeller_performance( ...
        performance_data, rpm_list, rpm, ...
        airspeed, prop_d, rho_reference);

    fprintf('\nPropeller at %.0f RPM:\n', rpm);
    fprintf('  Thrust: %.2f N\n', thrust);
    fprintf('  Shaft Power: %.2f W\n', shaft_power);
    fprintf('  Torque: %.3f Nm\n', torque);

end


%% 4. Test motor performance

voltage = 22.2;
rpm = 8000;

results = calculate_motor_performance(motor, rpm, voltage);

fprintf('\nMotor at %.0f RPM and %.1f V:\n', rpm, voltage);
fprintf('  Current: %.2f A\n', results.current);
fprintf('  Torque: %.3f Nm\n', results.torque);
fprintf('  Input Power: %.2f W\n', results.power_in);
fprintf('  Efficiency: %.2f%%\n', results.efficiency * 100);
fprintf('  Feasible: %d\n', results.feasible);


%% 5. Motor-propeller matching
% To be implemented next
%% Evaluate torque difference at two RPM values

rpm1 = 7000;
rpm2 = 8000;
delta_Q1 = match_motor_propeller(motor, rpm1, voltage, ...
    performance_data, rpm_list, airspeed, prop_d, rho_reference);
delta_Q2 = match_motor_propeller(motor, rpm2, voltage, ...
    performance_data, rpm_list, airspeed, prop_d, rho_reference);

fprintf("Value of Delta Q1: %d\n", delta_Q1);
fprintf("Value of Delta Q2: %d\n", delta_Q2);
%% 7. Search for motor-propeller torque equilibrium

rpm_test = 7000:500:11000;

rpm_prev = NaN;
delta_Q_prev = NaN;
rpm_lower_bound = NaN;
rpm_upper_bound = NaN;
for rpm = rpm_test
    delta_Q = match_motor_propeller(motor, rpm, voltage, ...
        performance_data, rpm_list, airspeed, prop_d, rho_reference);

    fprintf('RPM: %.0f | Delta Q: %.3f Nm\n', rpm, delta_Q);
  
    if (~isnan(delta_Q_prev) && ~isnan(delta_Q) ...
            && delta_Q_prev * delta_Q < 0) % detect sign change) 
        
        rpm_lower_bound = rpm_prev;
        rpm_upper_bound = rpm;
        break;
    end
    delta_Q_prev = delta_Q;
    rpm_prev = rpm;
end
fprintf('\nLower RPM bound: %.0f\n', rpm_lower_bound);
fprintf('Upper RPM bound: %.0f\n', rpm_upper_bound);
%% 8. Solve for motor-propeller equilibrium

% Check for valid bounds and find zero
if ~isnan(rpm_lower_bound) && ~isnan(rpm_upper_bound)

    torque_error = @(rpm) match_motor_propeller( ...
        motor, rpm, voltage, ...
        performance_data, rpm_list, airspeed, prop_d, rho_reference);
    rpm_equilibrium = fzero(torque_error, [rpm_lower_bound, rpm_upper_bound]);

    fprintf('Equilibrium RPM: %.2f\n', rpm_equilibrium);
else
    fprintf('No torque equilibrium bracket found.\n');
end

%% 9. Evaluate Operating Conditions
motor_results = calculate_motor_performance(motor, rpm_equilibrium, voltage);
[thrust_eq, shaft_power_eq, torque_prop_eq] = ...
    calculate_propeller_performance( ...
    performance_data, rpm_list, rpm_equilibrium, ...
    airspeed, prop_d, rho_reference);
fprintf("Equilibrium Current: %d\n", motor_results.current);
fprintf("Motor Torque: %d\n", motor_results.torque);
fprintf("Efficiency: %d\n", motor_results.efficiency);
fprintf("Feasibility: %d\n", motor_results.feasible);
fprintf("Torque Balance: %d\n", torque_prop_eq);
fprintf("Thrust @ Equilibrium: %d\n", thrust_eq);
fprintf("Shaft Power @ Equilibrium: %d\n", shaft_power_eq);
torque_residual = motor_results.torque - torque_prop_eq