# DBF 2026-2027 Optimizer 
Last Updated - 9/16/26

## List of variables accounted for by optimizer 
| Category | Parameter | Symbol | Units | Purpose |
| :--- | :--- | :--- | :--- | :--- |
| **Aerodynamics** | Wing Area | $S$ | $m^2$ | Calculates the total Lift and total Drag forces. |
| **Aerodynamics** | Aspect Ratio | $AR$ | $m$ / Unitless | Calculates induced drag. High Aspect Ratio drastically increases range. |
| **Aerodynamics** | Zero-Lift Drag | $C_{D0}$ | Unitless | The airframe's base friction. Crucial for calculating absolute top speed. |
| **Aerodynamics** | Oswald Efficiency | $e$ | Unitless | Efficiency factor (0.7-0.85) accounting for wingtip vortices and fuselage interference. |
| **Aerodynamics** | Max Lift Coefficient | $C_{L,max}$ | Unitless | Determines stall speed (the absolute minimum speed required to fly). |
| **Propulsion** | Motor Max Power | $P_{max}$ | $W$ | The maximum electrical Watts the motor can safely pull to fight drag at top speed. |
| **Propulsion** | Motor Efficiency | $\eta_{motor}$ | Unitless ($0.0 - 1.0$) | Percentage of electrical power successfully converted to physical shaft power. |
| **Propulsion** | Propeller Efficiency | $\eta_{prop}$ | Unitless ($0.0 - 1.0$) | Percentage of shaft rotation successfully converted into forward thrust. |
| **Propulsion** | Prop Pitch & Diameter | - | $m$ | The "gearing." Low pitch limits top speed; high pitch reduces low-speed thrust. |
| **Mass & Energy** | Total All-Up Weight | $W$ | $N$ | Total flying weight force (mass in $kg$ × $9.81$). High weight demands more lift, increasing drag. |
| **Mass & Energy** | Battery Capacity | $C$ | $mAh$ | Total milliampere-hours. When combined with voltage, calculates total flight time. |
| **Mass & Energy** | Battery Voltage | $V$ | $V$ | Dictates the motor's RPM ceiling and the total electrical power available. |
| **Environment** | Air Density | $\rho$ | $kg/m^3$ | Density of the air ($1.225 \text{ kg/m}^3$ at sea level). Affects drag, lift, and prop bite. |
| **Environment** | Gravity | $g$ | $m/s^2$ | Constant ($9.81 \text{ m/s}^2$) used to convert the plane's mass ($kg$) into downward force ($N$). |

## runmefirst.txt 
### This file lists all prerequisite python modules needed to run this program. 
Install the modules in this file by running them in terminal/bash!
