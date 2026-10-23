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