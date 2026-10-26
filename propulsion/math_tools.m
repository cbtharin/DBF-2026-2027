function tools = math_tools()
% Map human-readable names to your math functions
tools.rpmToRadSec = @convert_rpm_to_radsec;
tools.radSecToRpm = @convert_radsec_to_rpm;
tools.torque      = @calc_torque;
end

% --- Conversion Helpers ---
function rad_sec = convert_rpm_to_radsec(rpm)
rad_sec = rpm * (2 * pi / 60);
end

function rpm = convert_radsec_to_rpm(rad_sec)
rpm = rad_sec * (60 / (2 * pi));
end

% --- Equation Helpers ---
function T = calc_torque(power, omega)
% Calculates torque given shaft power (W) and angular velocity (rad/s)
if omega == 0
    T = 0;
else
    T = power / omega;
end
end
