function [rpm_list, static_data, prop_inches] = verify_static_data(filename)
%use OG file unless another is specified
if nargin < 1
    filename = fullfile("reference_data", "PER3_17x10E.dat");
end
% first we read the file
apc_text = fileread(filename);

% extracting the RPM section headings %
rpm_positions = strfind(apc_text, "PROP RPM")

%init array of all rpm values
rpm_list = zeros(1, length(rpm_positions));

% start 0 MPH PROP Extraction
for i=1:length(rpm_list)
    get_rpm = apc_text(rpm_positions(i): rpm_positions(i)+30);
    parts = split(get_rpm, "=");
    rpm_list(i) = str2double(parts(2));
end
rpm_list; %display select PROPS

% initialize the data array using length so that for any amount of rpm
% positions the parser will work
static_data = zeros(length(rpm_positions), 15);
size(static_data);


for i=1:length(rpm_positions) % processes each of the 14 sections
start_pos = rpm_positions(i);
if i < length(rpm_positions)
    end_pos = rpm_positions(i+1)-1;
else
    end_pos = length(apc_text); %last char of the file
end
% extracting row by row
section = apc_text(start_pos: end_pos);
section_line = splitlines(string(section)); % section line is an array of strings

% empty check example and finding static row to put in static matrix
    for j=1:length(section_line) % examines a single line (j) of text within a section
        row_data = sscanf(section_line(j), "%f");
        if ~isempty(row_data)
            if length(row_data) == 15 && row_data(1) == 0 % checks the result from sscanf() -> a col vector and validates length
                %disp("Line " + j + " of section " + i + " contains the static performance data");
                static_data(i, :) = row_data.'; % transpose the col to row vector
                break;
            end
        end
    
    end
end
static_data(1, :);
data_thrust = static_data(:, 11);

end
