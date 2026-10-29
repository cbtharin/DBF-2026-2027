%{
[rpm_list, static_data] = verify_static_data();
%out = sim(newest_singleprop); % runs simulink model automatically
%logs = out.logsout;

data_CT = static_data(:, 4);
data_CP = static_data(:, 5);
data_thrust = static_data(:, 11); % ensures that im comparing with data from the parser
rpm_list = rpm_list.'

%look into estimating wing loading
Ct = data_CT(4)
rpm_list(4)
static_data(4,11)
%rev_per_sec = rpm_list(4)./60;
rho;
prop_d
rev_per_sec = rpm_list./60;
rho_apc = data_thrust./(data_CT.*rev_per_sec.^2*prop_d^4) % approx 1.212
rho_reference = 1.212;
% pull coefficients


% comparing thrust
calculated_thrust = data_CT.* rev_per_sec.^2 * rho * prop_d^4
predicted_thrust = data_CT.* rev_per_sec.^2 .* rho_reference * prop_d^4
percent_error = abs(data_thrust() - predicted_thrust())./data_thrust() .* 100


% getting values from simulink
thrust_signal = logs.get("thrust");
thrust_values = thrust_signal.Values.Data;
sim_thrust = thrust_values(end)
rpm_signal = logs.get("RPM");
rpm_values = rpm_signal.Values.Data;
sim_rpm = rpm_values(end)
apc_thrust = interp1(rpm_list, data_thrust, sim_rpm)
thrust_adj = 100 * (rho/rho_reference)
thrust_error = abs(sim_thrust-apc_thrust)/apc_thrust *100

% interp1(rpm_list, data_thrust, sim_thrust)

%}
%% STATIC PROPELLER VERIFICATION
% Compare Simulink predictions against APC reference data.

%% 1. Initialize model parameters and run simulation

initialize_propulsion;

%% 2. Load APC static reference data

[rpm_list, static_data] = verify_static_data();

rpm_list = rpm_list(:);

data_CT = static_data(:, 4);
data_CP = static_data(:, 5);
data_thrust = static_data(:, 11);   % Newtons
data_torque = static_data(:, 10);   % N*m
data_power = static_data(:, 9);     % Watts

%% 3. Extract simulation results

thrust_signal = logs.get("thrust");
omega_signal = logs.get("angular_velocity");

sim_thrust = thrust_signal.Values.Data(end);
sim_omega = omega_signal.Values.Data(end);

% Convert rad/s to RPM
sim_rpm = sim_omega * 60 / (2*pi);

% Extract simulated torque and shaft power
torque_signal = logs.get("shaft_torque"); % gives simulink info but you need to extract a number
power_signal = logs.get("shaft_power");

sim_torque = torque_signal.Values.Data(end);
sim_power = power_signal.Values.Data(end);

%% 4. Interpolate APC reference data at simulated RPM

apc_thrust = interp1(rpm_list, data_thrust, sim_rpm);
apc_torque = interp1(rpm_list, data_torque, sim_rpm);
apc_power = interp1(rpm_list, data_power, sim_rpm);

%% 5. Adjust APC reference data for simulation air density

rho_reference = 1.212;  % Approximate density inferred from APC data
rho_sim = rho_reference; % Density used in the Simulink model
density_ratio = rho_sim / rho_reference;

apc_thrust_adjusted = apc_thrust * density_ratio;
apc_torque_adjusted = apc_torque * density_ratio;
apc_power_adjusted = apc_power * density_ratio;

%% 6. Calculate thrust error

thrust_error = ...
    abs(sim_thrust - apc_thrust_adjusted) ...
    / abs(apc_thrust_adjusted) * 100;

%% Torque and power verification

torque_error = ...
    abs(sim_torque - apc_torque_adjusted) ...
    / abs(apc_torque_adjusted) * 100;

power_error = ...
    abs(sim_power - apc_power_adjusted) ...
    / abs(apc_power_adjusted) * 100;
%% 7. Display thrust results

fprintf("\n====================================\n");
fprintf("   STATIC PROPELLER VERIFICATION\n");
fprintf("====================================\n");

fprintf("Simulated RPM:         %.1f RPM\n", sim_rpm);
fprintf("Simulation density:    %.3f kg/m^3\n", rho_sim);
fprintf("APC reference density: %.3f kg/m^3\n", rho_reference);

fprintf("\nTHRUST COMPARISON\n");
fprintf("Simulated thrust:      %.3f N\n", sim_thrust);
fprintf("APC reference thrust:  %.3f N\n", apc_thrust);
fprintf("APC at sim density:    %.3f N\n", apc_thrust_adjusted);

fprintf("Thrust error:          %.2f %%\n", thrust_error);

fprintf("\nTORQUE COMPARISON\n");
fprintf("Simulated torque:      %.3f N*m\n", sim_torque);
fprintf("APC reference torque:  %.3f N*m\n", apc_torque_adjusted);
fprintf("Torque error:          %.2f %%\n", torque_error);

fprintf("\nSHAFT POWER COMPARISON\n");
fprintf("Simulated power:       %.3f W\n", sim_power);
fprintf("APC reference power:   %.3f W\n", apc_power_adjusted);
fprintf("Power error:           %.2f %%\n", power_error);

fprintf("====================================\n");
n = rpm_list(:) / 60;

rho_apc = data_thrust ./ ...
    (data_CT .* n.^2 .* prop_d^4);

disp(rho_apc)
fprintf("Mean implied APC density: %.4f kg/m^3\n", mean(rho_apc));

% Linear Interpolation of Advance Ratios

%f = (J - J1)/(J2 - J1);
%Ct = Ct1 + f(Ct2 - Ct1);
% Using Interp1 function in MATLAB
%y_target = interp1(x_data, y_data, x_target);

