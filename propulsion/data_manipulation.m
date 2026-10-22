[rpm_list, static_data] = verify_static_data();

data_CT = static_data(:, 4);
data_thrust = static_data(:, 11); % ensures that im comparing with data from the parser
rpm_list = rpm_list.'

%look into estimating wing loading
data_CT(4)
rpm_list(4)
static_data(4,11)
%rev_per_sec = rpm_list(4)./60;
rho;
prop_d
rev_per_sec = rpm_list./60;
rho_apc = data_thrust./(data_CT.*rev_per_sec.^2*prop_d^4) % approx 1.212
rho_reference = 1.212;
% comparing thrust
calculated_thrust = data_CT.* rev_per_sec.^2 * rho * prop_d^4
predicted_thrust = data_CT.* rev_per_sec.^2 .* rho_reference * prop_d^4
percent_error = abs(data_thrust() - predicted_thrust())./data_thrust() .* 100
