rpm = [1000; 2000; 3000; 4000];
thrust = [1.083; 4.348; 9.816; 17.521];
%plot(rpm, thrust);
%xlabel("Horizonal");
%ylabel("Vertical");
apc_text = fileread("reference_data/PER3_17x10E.dat");
does_contain = contains(apc_text, "PROP RPM") %returns logic
check = rpm_positions(1)
extract_v = contains(lines(5), "0.00")
for i = 1:12
    fprintf("Line %d: %s\n", i, lines(i));
end

%getting thrust from each section_line
%{
static_row = sscanf(section_line(5), "%f"); %assumes the data always starts on line 5 so this is not the best practice
rpm_thrust(i) = static_row(11);
rpm_CT(i) = static_row(4);
rpm_CP(i) = static_row(5);
%}

rpm_thrust;
rpm_CP; 
rpm_CT;


% parser build with examples 

rpm = [1000; 2000; 3000; 4000];
thrust = [1.083; 4.348; 9.816; 17.521];
%plot(rpm, thrust);
%xlabel("Horizonal");
%ylabel("Vertical");
apc_text = fileread("reference_data/PER3_17x10E.dat");

% extracting the RPM section headings %
rpm_positions = strfind(apc_text, "PROP RPM")

%init arrays
rpm_list = zeros(1, length(rpm_positions));
rpm_thrust = zeros(1, length(rpm_positions));
rpm_CT = zeros(1, length(rpm_positions));
rpm_CP = zeros(1, length(rpm_positions));

% start 0 MPH PROP Extraction
for i=1:length(rpm_list)
    get_rpm = apc_text(rpm_positions(i): rpm_positions(i)+30);
    extract_rpm = strtrim(get_rpm);
    parts = split(get_rpm, "=");
    rpm_list(i) = str2double(parts(2));
end
rpm_list %display select PROPS

% trying to get a row of data %
first_section = apc_text(rpm_positions(1): rpm_positions(2)-1);
%subtract 1 bc the second postion is the first char of the next rpm heading
lines = splitlines(string(first_section));
static_row = sscanf(lines(5), "%f");
%static_thrust = static_row(11)

% example code line below
static_data = zeros(14, 15);
size(static_data);
data_thrust = static_data(:, 11);


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
                disp("Line " + j + " of section " + i + " contains the static performance data");
                static_data(i, :) = row_data.'; % transpose the col to row vector
                break;
            end
        end

    end
    %getting thrust from each section_line
    %{
static_row = sscanf(section_line(5), "%f"); %assumes the data always starts on line 5 so this is not the best practice
rpm_thrust(i) = static_row(11);
rpm_CT(i) = static_row(4);
rpm_CP(i) = static_row(5);
    %}
end
static_data(1, :);
data_thrust = static_data(:, 11);



%{
rpm_thrust;
rpm_CP; 
rpm_CT;
%}


% trying to get a row of data hardcoding %
first_section = apc_text(rpm_positions(1): rpm_positions(2)-1);
%subtract 1 bc the second postion is the first char of the next rpm heading
lines = splitlines(string(first_section));
static_row = sscanf(lines(5), "%f");
%static_thrust = static_row(11)