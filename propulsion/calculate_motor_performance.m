function results = calculate_motor_performance(motor, rpm, voltage)

% calculate motor torque constant (Nm/A)
Kt = 60/(2 * pi * motor.Kv); 
results.Kt = Kt;

%calculate back emf in (V)
back_emf = rpm/motor.Kv;
results.back_emf = back_emf;

% calculate motor current (A)
current = (voltage - back_emf) ./ motor.Rm; % eventually we will do multiple rpms simultaneously
results.current = current;

% calculate motor shaft torque (Nm)
torque = Kt .* (current - motor.I0);
results.torque = torque;

% calculate electrical power (W)
power_in = voltage .* current;
results.power_in = power_in;

%convert rpm to angular velocity (rad/s)
omega = rpm .* (2*pi/60);

% calculate mechanical power (W)
power_out = torque .* omega;
results.power_out = power_out;

% calculate motor efficiency
efficiency = power_out ./ power_in;
results.efficiency = efficiency;

% Check motor operating constraints
results.current_ok = current >= 0 & current <= motor.Imax;

results.power_ok = power_in >= 0 & power_in <= motor.Pmax;

% Check whether operating point is feasible
results.feasible = results.current_ok & results.power_ok;


end

% How the solver works
% 
% Suppose we select:
% 
% * Motor A
% * APC 17×10E propeller
% * A particular battery voltage
% * Airspeed of 20 m/s
% 
% MATLAB would:
% 
% 1. Estimate motor torque at a trial RPM.
% 2. Calculate propeller torque at that RPM.
% 3. Compare the two.
% 4. Adjust RPM until the torques balance.
% 5. Calculate thrust, current, and power at the resulting RPM.
