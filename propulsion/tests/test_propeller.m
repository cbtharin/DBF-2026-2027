function tests = test_propeller
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
[rpm_list, perf] = extract_rpm_performance('PER3_17x10E.dat');
testCase.TestData.rpm_list = rpm_list;
testCase.TestData.perf = perf;
testCase.TestData.D = 17 * 0.0254;
testCase.TestData.rho = 1.212;
end

function testKnownPoint7000(testCase)
d = testCase.TestData;
[T, P, Q] = calculate_propeller_performance(d.perf, d.rpm_list, 7000, 34, d.D, d.rho);
verifyEqual(testCase, T, 6.3802, 'RelTol', 1e-3);
verifyEqual(testCase, P, 351.7459, 'RelTol', 1e-3);
verifyEqual(testCase, Q, 0.4798, 'RelTol', 1e-3);
end

function testTorqueTimesOmegaIsPower(testCase)
d = testCase.TestData;
[~, P, Q] = calculate_propeller_performance(d.perf, d.rpm_list, 7000, 20, d.D, d.rho);
verifyEqual(testCase, Q * 2*pi*7000/60, P, 'RelTol', 1e-12);
end

function testThrustFallsWithAirspeed(testCase)
d = testCase.TestData;
T = arrayfun(@(v) calculate_propeller_performance(d.perf, d.rpm_list, 7000, v, d.D, d.rho), [10 20 30]);
verifyTrue(testCase, all(diff(T) < 0));
end

function testEfficiencyBelowOne(testCase)
d = testCase.TestData;
[T, P] = calculate_propeller_performance(d.perf, d.rpm_list, 7000, 20, d.D, d.rho);
verifyLessThan(testCase, T*20/P, 1);
end

function testOutOfRangeJ(testCase)
d = testCase.TestData;
verifyError(testCase, @() calculate_propeller_performance(d.perf, d.rpm_list, 7000, 60, d.D, d.rho), ...
    'prop:JOutOfRange');
end

% Next three functions protect a layer of the motor work

function testParserFindsAllSections(testCase)
d = testCase.TestData;
verifyEqual(testCase, d.rpm_list, 1000:1000:14000);
verifyEqual(testCase, numel(d.perf), 14);
verifyEqual(testCase, size(d.perf{d.rpm_list == 8000}), [30 15]);
end

function testExtractDimensions(testCase)
[D_in, pitch_in] = extract_dimensions('PER3_17x10E.dat');
verifyEqual(testCase, D_in, 17);
verifyEqual(testCase, pitch_in, 10);
end

function testSweepShapeAndColumns(testCase)
d = testCase.TestData;
R = propulsion_analysis(d.perf, d.rpm_list, [7000 8000], [10 20 30], d.D, d.rho);
verifyEqual(testCase, size(R), [6 5]);
verifyEqual(testCase, R(:,1).', [7000 7000 7000 8000 8000 8000]);
verifyEqual(testCase, R(:,2).', [10 20 30 10 20 30]);
[T, ~, ~] = calculate_propeller_performance(d.perf, d.rpm_list, 8000, 20, d.D, d.rho);
verifyEqual(testCase, R(5,3), T, 'RelTol', 1e-12);
end