function [rpm_list, performance_data] = extract_rpm_performance(filename)
%use OG file unless another is specified
if nargin < 1
    filename = fullfile("reference_data", "PER3_17x10E.dat");
end
% first we read the file
apc_text = fileread(filename);

% extracting the RPM section headings %
rpm_positions = strfind(apc_text, "PROP RPM");

%another way to do it positions = regexp(apc_text, 'PROP RPM\s*=\s*\d+', 'start');

%init array of all rpm values
rpm_list = zeros(1, length(rpm_positions));
performance_data = cell(length(rpm_positions), 1);
% start 0 MPH PROP Extraction
for i=1:length(rpm_list)
    get_rpm = apc_text(rpm_positions(i): rpm_positions(i)+30);
    parts = split(get_rpm, "=");
    rpm_list(i) = str2double(parts(2));
end

% lets separate into sections by rpm value
for i=1:length(rpm_positions) % processes each of the 14 sections
    section_data = zeros(0, 15);
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
        
        if length(row_data) == 15
            section_data(end+1, :) = row_data.';
        end  

    end
    performance_data{i} = section_data;
end
end

