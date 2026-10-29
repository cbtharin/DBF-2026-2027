function tests = test_motor_prop_match
tests = functiontests(localfunctions);
end

function setupOnce(testCase)

% Define reference motor
motor.Kv = 520;
motor.Rm = 0.016;
motor.I0 = 1.40;
motor.Imax = 100;
motor.Pmax = 2000;

% Load propeller data
% Find the propulsion folder relative to this test file
test_dir = fileparts(mfilename('fullpath'));
propulsion_dir = fileparts(test_dir);


% Build an absolute path to the APC data file
prop_file = fullfile(propulsion_dir, 'reference_data', 'PER3_17x10E.dat');
[rpm_list, perf] = extract_rpm_performance(prop_file);

% Solve reference configuration
results = solve_motor_prop_match( ...
    motor, 22.2, perf, rpm_list, 20, ...
    17 * 0.0254, 1.212);

% Save results for the individual tests
testCase.TestData.results = results;

end

function testTorqueEquilibrium(testCase)

results = testCase.TestData.results;

verifyLessThan(testCase, ...
    abs(results.motor_torque - results.prop_torque), 1e-3);

end

function testKnownEquilibriumRPM(testCase)

results = testCase.TestData.results;

verifyEqual(testCase, results.rpm, 10138.25, 'AbsTol', 1);

end

function testInfeasibleOperatingPoint(testCase)

results = testCase.TestData.results;

verifyFalse(testCase, results.feasible);

end