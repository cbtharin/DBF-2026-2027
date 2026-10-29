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