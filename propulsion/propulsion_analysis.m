function [analysis_results] = propulsion_analysis(performance_data, rpm_list, rpm_values, airspeed_values, diameter, rho_reference)
    
analysis_results = NaN(length(rpm_values) * length(airspeed_values), 5);

row = 1;
for i=1:length(rpm_values)
    for j=1:length(airspeed_values)
    % analysis_results(row, 1) = rpm_values(i);
    % airspeed = airspeed_values(j);
    % [thrust, shaft_power, torque] = calculate_propeller_performance(performance_data, rpm_list, rpm_values(i), airspeed, diameter, rho_reference);
    % analysis_results(row, 2) = airspeed; % airspeed
    % analysis_results(row, 3) = thrust; % thrust
    % analysis_results(row, 4) = shaft_power; % power
    % analysis_results(row, 5) = torque; % torque
    % row = row+1;
    analysis_results(row, 1) = rpm_values(i);

    airspeed = airspeed_values(j);
    analysis_results(row, 2) = airspeed;

    % Calculate advance ratio for this operating point
    rpm_index = find(rpm_list == rpm_values(i), 1);
    prop_table = performance_data{rpm_index};

    J_data = prop_table(:, 2);

    n = rpm_values(i) / 60;
    J_target = airspeed / (n * diameter);

    % Only calculate if the APC data supports this advance ratio
    if J_target >= min(J_data) && J_target <= max(J_data)

        [thrust, shaft_power, torque] = calculate_propeller_performance( ...
            performance_data, rpm_list, rpm_values(i), ...
            airspeed, diameter, rho_reference);

        analysis_results(row, 3) = thrust;
        analysis_results(row, 4) = shaft_power;
        analysis_results(row, 5) = torque;

    end

    row = row + 1;
    end
    
end

end


