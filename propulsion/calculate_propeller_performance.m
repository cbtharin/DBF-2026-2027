function [thrust, shaft_power, torque, APC_reference] = calculate_propeller_performance(performance_data, ...
    rpm_list, rpm, airspeed, diameter, rho_reference)

rpm_index = find(rpm_list == rpm); 
if isempty(rpm_index)
    error('Requested RPM value is not available in rpm_list.');
end
prop_table = performance_data{rpm_index};
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

% Interpolate coefficients and values (btwn advance ratios within a table)
Ct_target = interp1(J_data, Ct_data, J_target, 'linear');
Cp_target = interp1(J_data, Cp_data, J_target, 'linear');

thrust = Ct_target * rho_reference * n^2 * diameter^4;
shaft_power = Cp_target * rho_reference *n^3 * diameter^5;
omega = (2*pi * n);
torque = shaft_power / omega; % Motor Torque P = Q * Ω

end