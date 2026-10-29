function tests = test_motor
tests = functiontests(localfunctions);
end

function testMotorOvercurrent(testCase)

motor.Kv = 520;
motor.Rm = 0.016;
motor.I0 = 1.40;
motor.Imax = 100;
motor.Pmax = 2000;

results = calculate_motor_performance(motor, 10138.25, 22.2);

verifyGreaterThan(testCase, results.current, motor.Imax);
verifyFalse(testCase, results.feasible);

end