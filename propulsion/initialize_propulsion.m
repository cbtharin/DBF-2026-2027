prop_d = 0.559; % propeller diameter in meters
rho = 1.075; %*1/2 is pressure coeffiecient
kt = 0.11;
kp = 0.055;
mass = 10.8862;
g = 9.81;
rad = 2*pi;
%CQ = 0;
%CT =0;
%RPM = 0;


logs = out.logsout; %can print all signals as matrix value

fprintf("Number of logged signals: %d\n\n", logs.numElements);

for k = 1:logs.numElements
    signal = logs.getElement(k);

    fprintf("Signal %d\n", k)
    fprintf("  Name: %s\n", char(signal.Name))

    % Extract the path from the BlockPath object
    blockPath = signal.BlockPath.getBlock(1);
    fprintf("  Source: %s\n\n", blockPath)
end

%wing surface area = 
%coefficient of lift 
%coefficient of drag in unitless and a constant = 
%area = 
%20 to 30 seconds to get to max speed