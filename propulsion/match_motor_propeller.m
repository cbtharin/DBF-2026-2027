function delta_Q = match_motor_propeller(motor, rpm, voltage, ...
    performance_data, rpm_list, airspeed, prop_d, rho_reference)

 % OUTDATED THIS IS FULLY FINISHED INSIDE SOLVE_MOTOR_PROP_MATCH


% calculate motor performance
motor_results = calculate_motor_performance(motor, rpm, voltage);

%calculate prop performance
[thrust, shaft_power, prop_torque] = calculate_propeller_performance( ...
    performance_data, rpm_list, rpm, airspeed, prop_d, rho_reference);

Qm = motor_results.torque;
Qp = prop_torque;
delta_Q = Qm - Qp;

% we need to find rpm where delta_Q equals 0
% delta_Q is a function of RPM

% torque_equilibrium = fzero(delta_Q);
end