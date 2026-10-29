motor.name = "Scorpion SII-4025-520";

motor.Kv = 520;       % RPM/V
motor.Rm = 0.016;     % Ohms
motor.I0 = 1.40;      % Amps, measured at 10 V

motor.Imax = 100;     % Amps
motor.Pmax = 2000;    % Watts
motor.mass = 0.353;   % kg

rpm = 8000;
voltage = 22.2;

results = calculate_motor_performance(motor, rpm, voltage);
disp(results)