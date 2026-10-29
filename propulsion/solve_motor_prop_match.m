function results = solve_motor_prop_match(motor, voltage, performance_data, ...
    rpm_list, airspeed, prop_d, rho_reference)

    % Estimate maximum motor RPM i.e. no load RPM
    rpm_max_motor = motor.Kv * voltage;
    
    % Generate RPM values for torque equilibrium 
    % linspace generates 100 evenly spaced values from rpm_min to rpm_max
    % allows the algorithm to adapt to different motors
    % need min and max check to ensure apc has data for requested rpm
    rpm_min = max(1000, min(rpm_list));
    rpm_max = min(rpm_max_motor, max(rpm_list));
    rpm_test = linspace(rpm_min, rpm_max, 100);
    
    if(rpm_min >= rpm_max)
    results.rpm = NaN;
    results.feasible = false;
    return;
    end
    %% 7. Search for motor-propeller torque equilibrium

    rpm_prev = NaN;
    delta_Q_prev = NaN;
    rpm_lower_bound = NaN;
    rpm_upper_bound = NaN;
    for rpm = rpm_test
        delta_Q = match_motor_propeller(motor, rpm, voltage, ...
            performance_data, rpm_list, airspeed, prop_d, rho_reference);

        fprintf('RPM: %.0f | Delta Q: %.3f Nm\n', rpm, delta_Q);

        if (~isnan(delta_Q_prev) && ~isnan(delta_Q) && ...
        delta_Q_prev * delta_Q < 0) % detect sign change) 

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
        results.rpm = rpm_equilibrium;
    else
        results.feasible = false;
        results.rpm = Nan;
        fprintf('No torque equilibrium bracket found.\n');
        return;
    end
    motor_results = calculate_motor_performance(motor, rpm_equilibrium, voltage);
    [thrust, shaft_power, prop_torque] = calculate_propeller_performance(performance_data, rpm_list, ...
        rpm_equilibrium, airspeed, prop_d, rho_reference);

    % Assigning fields
    % from motor
    results.current = motor_results.current;
    results.motor_torque = motor_results.torque;
    results.efficiency = motor_results.efficiency;
    results.feasible = motor_results.feasible;

    % from prop
    results.shaft_power = shaft_power;
    results.prop_torque = prop_torque;
    results.thrust = thrust;
end

 

