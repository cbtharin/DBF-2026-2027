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