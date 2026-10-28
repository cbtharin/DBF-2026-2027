function [prop_inches, pitch_inches] = extract_dimensions(filename)

% extracting prop dimensions 

result = regexp(filename, '\d+(\.\d+)?x\d+(\.\d+)?', 'match', 'once');
if isempty(result)
    error("Invalid propeller filename.");
end
dimensions = split(result, "x");
dimensions = str2double(dimensions);

% Normalize APC abbreviated dimensions
for i = 1:2
    value = dimensions(i);

    if value >= 3000
        dimensions(i) = value / 1000;

    elseif value >= 300
        dimensions(i) = value / 100;

    elseif value >= 100
        dimensions(i) = value / 10;

    elseif value >= 25
        dimensions(i) = value / 10;
    end
end

prop_inches = dimensions(1);
prop_inches = prop_inches * 0.0254;
pitch_inches = dimensions(2);

end