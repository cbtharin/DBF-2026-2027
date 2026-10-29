function [thrust, shaft_power, torque] = calculate_propeller_performance(performance_data, ...
    rpm_list, rpm, airspeed, diameter, rho_reference)


% Find RPM sections surrounding requested RPM
rpm_below = rpm_list(rpm_list <= rpm);
rpm_above = rpm_list(rpm_list >= rpm);

% Check whether requested RPM is within available data
if isempty(rpm_below) || isempty(rpm_above)
    error('Requested RPM is outside the APC data range.');
end

rpm_lower = max(rpm_below);
rpm_upper = min(rpm_above);

% find returns the position where to condition is true
rpm_index_l = find(rpm_list == rpm_lower);
rpm_index_u = find(rpm_list == rpm_upper);
prop_table_upper = performance_data{rpm_index_u};
prop_table_lower = performance_data{rpm_index_l};

% Extract the upper/lower coefficients
J_data_lower = prop_table_lower(:, 2);
Ct_data_lower = prop_table_lower(:, 4);
Cp_data_lower = prop_table_lower(:, 5);
thrust_data_lower = prop_table_lower(:, 11);
shaft_power_data_lower = prop_table_lower(:, 9);
torque_data_lower = prop_table_lower(:, 10);

J_data_upper = prop_table_upper(:, 2);
Ct_data_upper = prop_table_upper(:, 4);
Cp_data_upper = prop_table_upper(:, 5);
thrust_data_upper = prop_table_upper(:, 11);
shaft_power_data_upper = prop_table_upper(:, 9);
torque_data_upper = prop_table_upper(:, 10);


% Extract the APC coefficients
J_data = prop_table(:, 2);
Ct_data = prop_table(:, 4);
Cp_data = prop_table(:, 5);
thrust_data = prop_table(:, 11);
shaft_power_data = prop_table(:, 9);
torque_data = prop_table(:, 10);

% Convert RPM to revolutions per second
n = rpm/60;
n_lower = rpm_lower/ 60;
n_upper = rpm_upper/60;

% Calculate advance ratio and range check 
J_target_lower = airspeed / (n_lower * diameter);
J_target_upper = airspeed / (n_upper * diameter);
if J_target_lower < min(J_data_lower) || J_target_lower > max(J_data_lower) ...
        || J_target_upper < min(J_data_upper) || J_target_upper > max(J_data_upper)
    error('Requested advance ratio is outside the APC data range.');
end

% Interpolate coefficients and values (btwn advance ratios within a table)
% J is independent var. First step to two stage interpolation 
Ct_target_lower = interp1(J_data_lower, Ct_data_lower, J_target_lower, 'linear');
Ct_target_upper = interp1(J_data_upper, Ct_data_upper, J_target_upper, 'linear');
Cp_target_lower = interp1(J_data_lower, Cp_data_lower, J_target_lower, 'linear');
Cp_target_upper = interp1(J_data_upper, Cp_data_upper, J_target_upper, 'linear');

% rpm is now independent var. Second part of Two Stage Interpolation
if (~(rpm_lower == rpm_upper))
Ct_target = interp1([rpm_lower, rpm_upper], [Ct_target_lower, Ct_target_upper], rpm, 'linear' );
Cp_target = interp1([rpm_lower, rpm_upper], [Cp_target_lower, Cp_target_upper], rpm, 'linear' );

else %no interpolation necessary
    Ct_target = Ct_target_lower;
    Cp_target = Cp_target_lower;
end

% final calcs use requested rpm not upper/lower
thrust = Ct_target * rho_reference * n^2 * diameter^4;
shaft_power = Cp_target * rho_reference *n^3 * diameter^5;
omega = (2*pi * n);
torque = shaft_power / omega; % Motor Torque P = Q * Ω

end