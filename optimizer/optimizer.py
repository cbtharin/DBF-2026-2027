import random
import math
import matplotlib.pyplot as plt

# 1. Setup lists to save the data of planes that successfully fly
speeds = []
ranges = []

print("Simulating 50,000 random RC planes... please wait.")

# 2. Start a giant loop to guess 50,000 different planes
for i in range(50000):
    
    # --- A. GUESS RANDOM VARIABLES ---
    # Pick a random number between a minimum and maximum limit
    S = random.uniform(0.2, 0.8)             # Wing Area (m^2)
    AR = random.uniform(4.0, 15.0)           # Aspect Ratio
    C_batt_Ah = random.uniform(2.0, 10.0)    # Battery Capacity (Ah)
    P_motor_max = random.uniform(200, 1500)  # Motor Power (Watts)

    # --- B. CONSTANTS ---
    rho = 1.225         # Air density
    g = 9.81            # Gravity
    Cd0 = 0.02          # Base drag
    e = 0.85            # Wing efficiency
    Cl_max = 1.2        # Max lift before stalling
    eta_motor = 0.80    # Motor efficiency
    eta_prop = 0.70     # Prop efficiency
    V_batt = 14.8       # Battery Voltage

    # --- C. CALCULATE WEIGHT ---
    m_wing = (0.5 * S) + (0.02 * (AR ** 1.5))
    m_batt = (V_batt / 3.7) * C_batt_Ah * 0.025
    m_motor = P_motor_max / 3000.0
    m_empty = 1.0
    
    W = (m_empty + m_wing + m_batt + m_motor) * g # Total weight in Newtons

    # --- D. CHECK CONSTRAINTS (Does it crash?) ---
    span = math.sqrt(S * AR)
    V_stall = math.sqrt((2 * W) / (rho * S * Cl_max))
    
    # If the wingspan is over 1.83m (6ft) OR stall speed is too fast, skip this plane
    if span > 1.83 or V_stall > 15.0:
        continue # 'continue' skips the rest of the loop and starts the next guess

    # --- E. CALCULATE SPEED AND RANGE ---
    # Max Speed
    P_shaft = P_motor_max * eta_motor * eta_prop
    V_max = ( (2 * P_shaft) / (rho * S * Cd0) ) ** (1/3)
    
    # Max Range (Cruising Efficiency)
    Cl_cruise = math.sqrt(Cd0 * math.pi * e * AR)
    Cd_cruise = 2 * Cd0
    V_cruise = math.sqrt((2 * W) / (rho * S * Cl_cruise))
    
    Drag_cruise = W / (Cl_cruise / Cd_cruise)
    P_elec_cruise = (Drag_cruise * V_cruise) / (eta_motor * eta_prop)
    
    Flight_time_seconds = (V_batt * C_batt_Ah * 3600) / P_elec_cruise
    Range_meters = V_cruise * Flight_time_seconds

    # --- F. SAVE THE WINNERS ---
    speeds.append(V_max)
    ranges.append(Range_meters / 1000) # Divide by 1000 to save as Kilometers

# 3. Graph the results
plt.figure(figsize=(8, 6))
plt.scatter(speeds, ranges, c="blue", alpha=0.3, s=5) # s=5 makes the dots small
plt.title("Monte Carlo Optimizer: 50,000 Random Planes")
plt.xlabel("Max Speed (m/s)")
plt.ylabel("Max Range (Kilometers)")
plt.grid(True)
plt.show()