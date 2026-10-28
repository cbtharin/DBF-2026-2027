%math = math_tools();

% Load APC performance data
[rpm_list, performance_data] = extract_rpm_performance();

% Define the operating condition
airspeed = 34;
airspeed_mph = [0, 10, 20, 30, 35, 40, 50, 60, 70, 76.1, 80];
airspeed_values = airspeed_mph * 0.44704;
rpm = 7000;
diameter = 17 * 0.0254;% Test propeller coefficient interpolation

% generate performance predictions
[analysis_results] = propulsion_analysis(performance_data, rpm_list, rpm_values, airspeed_values, diameter, rho_reference);

% Locate the requested RPM this case is 7000
rpm_index = find(rpm_list == rpm); 
%inside gives a logical array and find returns the location where array = 1

% Retrieve its performance table
prop_table = performance_data{rpm_index};
% init validation matrix
validation_errors = zeros(size(analysis_results, 1), 3); % the 1 represents which dimension to measure

% Extract the APC coefficients
J_data = prop_table(:, 2);
Ct_data = prop_table(:, 4);
Cp_data = prop_table(:, 5);
thrust_data = prop_table(:, 11);
shaft_power_data = prop_table(:, 9);
torque_data = prop_table(:, 10);

% Convert RPM to revolutions per second
n = rpm / 60;

% Calculate advance ratio and range check 
J_target = airspeed / (n * diameter);
if J_target < min(J_data) || J_target > max(J_data)
    error('Requested advance ratio is outside the APC data range.');
end

% Interpolate coefficients and values (values that are not between
% thousands)
Ct_target = interp1(J_data, Ct_data, J_target, 'linear');
Cp_target = interp1(J_data, Cp_data, J_target, 'linear');
APC_thrust = interp1(J_data, thrust_data, J_target, 'linear');
APC_power = interp1(J_data, shaft_power_data, J_target, 'linear');
APC_torque = interp1(J_data, torque_data, J_target, 'linear');
APC_reference = [apc_thrust, apc_power, apc_torque];

fprintf('Advance ratio: %.4f\n', J_target);
fprintf('Thrust coefficient: %.5f\n', Ct_target);
fprintf('Power coefficient: %.5f\n', Cp_target);

%fprintf('Available J range: %.4f to %.4f\n', J_min, J_max);
fprintf('Requested J: %.4f\n', J_target);

[thrust, shaft_power, torque] = calculate_propeller_performance( ...
    performance_data, rpm_list, rpm, airspeed, diameter, rho_reference);
fprintf("Thrust: %.3f\n", thrust);
fprintf("Shaft Power: %.3f\n", shaft_power);
fprintf("Motor Torque: %.3f N*m\n", torque);

% Error Calc
thrust_error = ...
    abs(thrust - APC_thrust) / abs(APC_thrust) * 100;
fprintf("Thrust Error: %.2f%%\n", thrust_error);
s_power_error = ...
    abs(shaft_power - APC_power) / abs(APC_power) * 100;
fprintf("Shaft Power Error: %.2f%%\n", s_power_error);
torque_error = abs(torque - APC_torque) / abs(APC_torque) * 100;
fprintf("Torque Error: %.2f%%\n", torque_error);

% looking for validation errors for thrust, power, and torque
for i = 1: size(analysis_results, 1) % explicitly looks for num of rows
    rpm = analysis_results(i, 1);
    airspeed = analysis_results(i, 2);
    thrust = analysis_results(i, 3);
    shaft_power = analysis_results(i, 4);
    torque = analysis_results(i, 5);

    % Find corresponding rpm and retrieve its performance table
    rpm_index = find(rpm_list == rpm);
    prop_table = performance_data{rpm_index};

    % Extract the APC coefficients
    J_data = prop_table(:, 2);
    thrust_data = prop_table(:, 11);
    shaft_power_data = prop_table(:, 9);
    torque_data = prop_table(:, 10);

    % Convert RPM to revolutions per second
    n = rpm / 60;

    % Calculate advance ratio and range check 
    J_target = airspeed / (n * diameter);
    if J_target < min(J_data) || J_target > max(J_data)
        error('Requested advance ratio is outside the APC data range.');
    end

    % Interpolate coefficients and values (values that are not between
    % thousands)
    APC_thrust = interp1(J_data, thrust_data, J_target, 'linear');
    APC_power = interp1(J_data, shaft_power_data, J_target, 'linear');
    APC_torque = interp1(J_data, torque_data, J_target, 'linear');

    % Error Calc
    thrust_error = ...
        abs(thrust - APC_thrust) / abs(APC_thrust) * 100;
    s_power_error = ...
        abs(shaft_power - APC_power) / abs(APC_power) * 100;
    torque_error = abs(torque - APC_torque) / abs(APC_torque) * 100;

    validation_errors(i, 1) = thrust_error
    validation_errors(i, 2) = s_power_error
    validation_errors(i, 3) = torque_error
end

max_errors = max(validation_errors, [], 1);
mean_errors = mean(validation_errors, 1);

fprintf('\nAPC Validation Summary\n');
fprintf('Operating conditions tested: %d\n', size(validation_errors, 1));

fprintf('Maximum errors [Thrust, Power, Torque] (%%):\n');
disp(max_errors);

fprintf('Mean errors [Thrust, Power, Torque] (%%):\n');
disp(mean_errors);


