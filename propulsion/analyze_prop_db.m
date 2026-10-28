prop_files = dir(fullfile("data", "apc-performance-files", "PERFILES2", "PER3_*.dat"));
%prop_files = prop_files(1:min(3, length(prop_files)));
num_files = length(prop_files);

% Define operating conditions
airspeed_mph = [0, 10, 20, 30, 35, 40, 50, 60, 70, 76.1, 80];
airspeed_values = airspeed_mph * 0.44704;

rho_reference = 1.212;

% each prop will have its own entry containing the things below
prop_results = struct('filename', {}, 'diameter', {}, 'pitch', {}, 'performance', {});

for i = 1:length(prop_files)

    % Get filename and full file path
    filename = prop_files(i).name;
    filepath = fullfile(prop_files(i).folder, filename);

    % Extract propeller dimensions
    [diameter, pitch] = extract_dimensions(filename);
    fprintf('%s -> %.2f x %.2f inches\n', ...
        filename, prop_inches, pitch_inches);[prop_inches, pitch_inches] = extract_dimensions(filename);

    %diameter = prop_inches * 0.0254;

    % Extract available RPM data
    [rpm_list, performance_data] = extract_rpm_performance(filepath);

    % Run the performance sweep
    analysis_results = propulsion_analysis( ...
        performance_data, rpm_list, rpm_list, ...
        airspeed_values, diameter, rho_reference);

    % Store results
    prop_results(i).filename = filename;
    prop_results(i).diameter = prop_inches;
    prop_results(i).pitch = pitch_inches;
    prop_results(i).performance = analysis_results;
end

for i = 1:length(prop_results)
    fprintf('%s: %.3f x %.3f inches, %d operating conditions\n', ...
        prop_results(i).filename, ...
        prop_results(i).diameter, ...
        prop_results(i).pitch, ...
        size(prop_results(i).performance, 1));
end