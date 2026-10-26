function [prop_inches, pitch_inches] = extract_dimensions(filename)

% extracting prop dimensions 

result = regexp(filename, '\d+(\.\d+)?x\d+(\.\d+)?', 'match', 'once');
if isempty(result)
    error("Invalid propeller filename.");
end
dimensions = split(result, "x");
dimensions = str2double(dimensions);
prop_inches = dimensions(1);
pitch_inches = dimensions(2);

end