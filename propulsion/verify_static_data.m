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
lines = splitlines(string(first_section))
static_row = sscanf(lines(5), "%f")
%static_thrust = static_row(11)

for i=1:length(rpm_positions)
start_pos = rpm_positions(i);
if i < length(rpm_positions)
    end_pos = rpm_positions(i+1)-1;
else
    end_pos = length(apc_text); %last char of the file
end
% extracting row by row
section = apc_text(start_pos: end_pos);
section_line = splitlines(string(section));
%getting thrust from each section_line
static_row = sscanf(section_line(5), "%f");
rpm_thrust(i) = static_row(11);
rpm_CT(i) = static_row(4);
rpm_CP(i) = static_row(5);
end
rpm_thrust
rpm_CP
rpm_CT
