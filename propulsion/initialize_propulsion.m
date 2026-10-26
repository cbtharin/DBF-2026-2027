%prop_d = 0.559; % propeller diameter in meters
prop_inches = 17;
prop_d = .0254 * prop_inches; %converts in to m 
rho = 1.075; %*1/2 is pressure coefficient
% kt = 0.11;
% kp = 0.055;
mass = 10.8862;
g = 9.81;
rad = 2*pi;

desired_rpm;

% get file specs
prop_file = fullfile("reference_data", "PER3_17x10E.dat");
[prop_inches, pitch_inches] = extract_dimensions(prop_file);


[rpm_list, static_data] = verify_static_data(prop_file);

data_CT = static_data(:, 4);

data_CP = static_data(:, 5);

rho_reference = 1.212;

%apc data is what the manufacturer predicts its purpose if reference
%prediction
%uiuc is experimental data its purpose is real world validation
%our model's purpose is to evaluate our theoretical calculations to make
%meaningful comparisons to the data sets

out = sim(newest_singleprop); % runs simulink model automatically

logs = out.logsout; %can print all signals as matrix value

fprintf("Number of logged signals: %d\n\n", logs.numElements);

for k = 1:logs.numElements
    signal = logs.getElement(k);

    fprintf("Signal %d\n", k);
    fprintf("  Name: %s\n", char(signal.Name));

    % Extract the path from the BlockPath object
    blockPath = signal.BlockPath.getBlock(1);
    %fprintf("  Source: %s\n\n", blockPath);
end

%wing surface area = 
%coefficient of lift 
%coefficient of drag in unitless and a constant = 
%area = 
%20 to 30 seconds to get to max speed

% Chosen Motor: Scorpion SII-4025-520-KV
Kv = 520; %RPM/V
battery = 65; %LiPo
nominal_v = 22.2; %A
max_cont_current = 100; %A
rated_power = 2000; %W