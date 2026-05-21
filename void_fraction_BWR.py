import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# Save figures next to the script in a figures/ subfolder.
_FIG_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")
os.makedirs(_FIG_DIR, exist_ok=True)

def _save_fig(name):
    """Save current figure to figures/<name>.png at 150 DPI."""
    plt.savefig(os.path.join(_FIG_DIR, f"{name}.png"), dpi=150, bbox_inches="tight")

# --- 1. Reactor Parameters ---
P_sys = 7.0e6           # System Pressure (Pa) 
Q_thermal = 4500e6      # Total Core Power (W) 
N_rods = 60000          # Number of rods
L_rod = 4.0             # Rod length (m)
D_rod = 1.1e-2          # Rod outer clad diameter (m)
v_inlet = 9.0           # Inlet flow velocity (m/s) 
pitch = 1.26e-2         # Pitch (m) 
T_inlet = 275.0         # Inlet temperature (C) 
dz = 0.01               # Step size

# Geometry Calculations
A_flow_subchannel = (pitch**2) - (np.pi * (D_rod/2)**2) 
P_wetted = np.pi * D_rod 
D_hyd = 4 * A_flow_subchannel / P_wetted 
D_e = D_hyd 

# --- 2. Material Properties at 7.0 MPa ---
rho_f = 739.72     # Liquid density (kg/m3)
rho_g = 36.525     # Vapor density (kg/m3)
mu_f = 9.1266e-05   # Liquid viscosity (Pa.s)
mu_g = 1.8889e-05   # Vapor viscosity (Pa.s)
sigma = 0.017633    # Surface tension (N/m)
P_crit = 22.064e6 # Critical Pressure (Pa)
g_const = 9.81    # Gravity





# --- 3. Correlations (Chexal-Lellouche & Friedel) ---

# 3.1 Advanced Bundle Friction (Cheng & Todreas with Geometric Correction)
# Sources:
# Eq 9.106 (Bundle Aggregation)
# Eq 9.109a (Geometry Coefficients)
# Table 9.5b (Turbulent Coefficients)

def get_friction_factor_cheng_todreas_bundle(Re_bundle):
    # --- 1. Define Geometry Parameters ---
    P = 1.26   # Pitch (cm)
    D = 1.10   # Diameter (cm)
    P_D_ratio = P / D
    
    # Assumption: Gap to wall is exactly P/2 (Standard symmetric assumption)
    
    # --- 2. Calculate Hydraulic Diameters (D_e) for each type ---
    
    # Interior Subchannel (Standard unit cell)
    A_int = (P**2) - (np.pi * D**2)/4.0
    Pw_int = np.pi * D
    De_int = 4.0 * A_int / Pw_int
    
    # Edge Subchannel (Rectangular: P x P/2)
    # Area = Box - Half Rod
    A_edge = (P * (P/2.0)) - (np.pi * D**2)/8.0
    # Wetted Perimeter = Half Rod + Wall Length (P)
    Pw_edge = (np.pi * D)/2.0 + P
    De_edge = 4.0 * A_edge / Pw_edge
    
    # Corner Subchannel (Square: P/2 x P/2)
    # Area = Corner Box - Quarter Rod
    A_corner = ((P/2.0)**2) - (np.pi * D**2)/16.0
    # Wetted Perimeter = Quarter Rod + Wall Length (P/2 + P/2)
    Pw_corner = (np.pi * D)/4.0 + P
    De_corner = 4.0 * A_corner / Pw_corner

    # --- CALCULATE TRUE BUNDLE HYDRAULIC DIAMETER ---
    # Weighted sum for the whole 8x8 assembly
    N_int = 36.0
    N_edge = 24.0
    N_corner = 4.0
    
    # Total Flow Area of the Bundle
    A_bundle_total = (N_int * A_int) + (N_edge * A_edge) + (N_corner * A_corner)
    
    # Total Wetted Perimeter of the Bundle
    Pw_bundle_total = (N_int * Pw_int) + (N_edge * Pw_edge) + (N_corner * Pw_corner)
    
    # True Hydraulic Diameter of the Bundle
    De_bundle = 4.0 * A_bundle_total / Pw_bundle_total
    
    # ... (Continue with normalization D_prime_eb = De_bundle / D) ...
    
    # --- 3. Calculate C' coefficients for each type (Table 9.5b) ---
    val = P_D_ratio - 1.0
    
    # Interior Coefficients
    C_int = 0.1339 + 0.09059*val - 0.09926*(val**2)
    
    # Edge Coefficients
    C_edge = 0.1430 + 0.04199*val - 0.04428*(val**2)
    
    # Corner Coefficients
    C_corner = 0.1452 + 0.02681*val - 0.03411*(val**2)
    
    # --- 4. Aggregate using Equation 9.106 ---
    # C'_b = D'_eb * [ Sum ( S_i * (D'_ei / D'_eb)^(n/(2-n)) * (C'_fi / D'_ei)^(1/(n-2)) ) ]^(n-2)
    # D'_ei = De_i / D (Normalized diameter)
    # C'_fi = C_i (The coefficient we just calculated)
    
    n = 0.20 
    
    # Normalize Diameters by Rod Diameter D (as per textbook convention D')
    D_prime_eb = De_bundle / D
    D_prime_int = De_int / D
    D_prime_edge = De_edge / D
    D_prime_corner = De_corner / D
    
    # Subchannel Area Fractions (S_i) for 8x8 Assembly
    S_int = 36.0 / 64.0
    S_edge = 24.0 / 64.0
    S_corner = 4.0 / 64.0
    
    # Exponents for Eq 9.106
    exp1 = n / (2.0 - n)       # Exponent for Diameter Ratio
    exp2 = 1.0 / (n - 2.0)     # Exponent for C/D Ratio
    exp_final = n - 2.0        # Final Exponent
    
    # Summation Terms
    term_int = S_int * ((D_prime_int/D_prime_eb)**exp1) * ((C_int/D_prime_int)**exp2)
    term_edge = S_edge * ((D_prime_edge/D_prime_eb)**exp1) * ((C_edge/D_prime_edge)**exp2)
    term_corner = S_corner * ((D_prime_corner/D_prime_eb)**exp1) * ((C_corner/D_prime_corner)**exp2)
    
    sum_total = term_int + term_edge + term_corner
    
    # Final Bundle Coefficient C'_bL (or C'_bT for turbulent)
    C_bundle = D_prime_eb * (sum_total ** exp_final)
    
    # --- 5. Final Friction Factor ---
    # f = C'_bundle / Re^n
    f_bundle = C_bundle * (Re_bundle**-n)
    
    return f_bundle

# 3.2 Friedel Two-Phase Multiplier (Eq 11.99)
def calculate_friedel_multiplier(x, G, D_e):
    # Homogeneous Density
    rho_m = 1.0 / ( (x / rho_g) + ((1 - x) / rho_f) )
    
    # Froude and Weber Numbers
    Fr = (G**2) / (g_const * D_e * rho_m**2)
    We = (G**2 * D_e) / (sigma * rho_m)
    
    # Single-phase Friction Factors
    Re_fo = (G * D_e) / mu_f
    f_fo = get_friction_factor_cheng_todreas_bundle(Re_fo)
    
    Re_go = (G * D_e) / mu_g
    f_go = get_friction_factor_cheng_todreas_bundle(Re_go)
    
    # E, F, H Terms
    E = (1 - x)**2 + (x**2) * ( (rho_f * f_go) / (rho_g * f_fo) )
    F = (x**0.78) * ((1 - x)**0.224)
    H = (rho_f / rho_g)**0.91 * (mu_g / mu_f)**0.19 * (1 - (mu_g / mu_f))**0.7
    
    # Safety for tiny numbers (approaching single phase)
    if Fr < 1e-9 or We < 1e-9: return 1.0 
    
    phi_sq = E + (3.24 * F * H) / ( (Fr**0.0454) * (We**0.035) )
    return phi_sq

# 3.3 Chexal-Lellouche 
def calculate_chexal_lellouche(x_quality, G, P_sys):
    if x_quality <= 0:
        return 0.0
    if x_quality >= 0.99:
        return 1.0
    j_g = (x_quality * G) / rho_g
    j_f = ((1 - x_quality) * G) / rho_f
    j_total = j_g + j_f
    
    Re_g = (rho_g * j_g * D_e) / mu_g
    Re_f = (rho_f * j_f * D_e) / mu_f
    
    Re_g_val = Re_g if (Re_g > Re_f or Re_g < 0) else Re_f
    
    # Constants
    term_A1 = np.exp(max(-85, min(85, (-Re_g_val / 60000))))
    A1 = 1 / (1 + term_A1)
    B1 = min(0.8, A1)
    
    D2 = 0.09144
    C7 = (D2 / D_e)**0.6
    
    if C7 < 1.0: 
        C8 = C7 / (1 - C7)
    else:
        C8 = 1000 
        
    if C7 >= 1:
        C4 = 1.0
    else:
        C4 = 1 / (1 - np.exp(-C8))

    C5 = (150 * (rho_g / rho_f))**0.5
    C6 = C5 / (1 - C5) 
        
    density_ratio = rho_f / rho_g
    if density_ratio >= 18:
        if C5 >= 1:
            C2 = 1.0
        elif C6 >= 85:
            C2 = 1.0
        else:
            C2 = 1 / (1 - np.exp(-C6))
    else:
        term_log = np.log(max(1.00001, density_ratio))
        C2 = 0.4757 * (term_log)**0.7

    D1 = 0.0381
    term_c10_1 = 2 * np.exp((abs(Re_f)/350000)**0.4) - 1.75 * (abs(Re_f)**0.035) * np.exp((-abs(Re_f)/60000) * (D1/D_e)**2)
    term_c10_2 = (D1/D_e)**0.1 * abs(Re_f)**0.001
    C10 = term_c10_1 + term_c10_2
    
    B2 = 1 / ((1 + 0.05 * abs(Re_f) / 350000)**0.4)
    
    if j_f > 0 and j_g > 0:
        C3 = max(0.50, 2 * np.exp(-abs(Re_f)/300000))
    else:
        C3 = 2 * (C10 / 2)**B2 

    def residual(alpha):
        if alpha <= 1e-6: alpha = 1e-6
        if alpha >= 0.999: alpha = 0.999

        if Re_g >= 0:
            C1 = (1 - alpha)**B1
        else:
            C1 = (1 - alpha)**0.5
            
        C_p = abs((4 * P_crit**2) / (P_sys * (P_crit - P_sys)))
        
        if alpha * C_p < 170:
            L_n = 1 - np.exp(-alpha * C_p)
        else:
            L_n = 1.0
            
        if C_p < 170:
            L_d = 1 - np.exp(-C_p)
        else:
            L_d = 1.0
        
        L_factor = L_n / L_d if L_d > 1e-9 else 1.0
        
        K0 = B1 + (1 - B1) * (rho_g / rho_f)**0.25
        r = (1 + 1.57 * (rho_g / rho_f)) / (1 - B1)
        C0 = L_factor / (K0 + (1 - K0) * alpha**r)
        
        term_drift = ((rho_f - rho_g) * g_const * sigma / (rho_f**2))**0.25
        V_vj = 1.41 * term_drift * C1 * C2 * C3 * C4
        
        return alpha * (C0 * j_total + V_vj) - j_g

    alpha_guess = 1 / (1 + ((1-x_quality)/x_quality) * (rho_g/rho_f))
    solution = fsolve(residual, alpha_guess)
    return solution[0]

# --- 4. Main Calculation Loops (Thermodynamics + Void Fraction) ---
G = rho_f * v_inlet
Q_channel = Q_thermal / N_rods
m_dot = G * A_flow_subchannel 
h_inlet = 1210.4e3 
h_sat_f = 1267.7e3 
h_sat_g = 2772.6e3 
h_fg = h_sat_g - h_sat_f

z_axis = np.arange(0, L_rod + dz, dz)

# Dictionary to help automate the pressure loop later
# We will store (h, x, alpha) in these lists
cases_data = {
    'Constant': {'h': [], 'x': [], 'alpha': []},
    'Sinusoidal': {'h': [], 'x': [], 'alpha': []}
}

# --- A. Thermodynamics & Void Calculation Loop ---

# CASE 1: Constant Heat Flux
for z in z_axis:
    h_z = h_inlet + (Q_channel / m_dot) * (z / L_rod)
    x_eq = (h_z - h_sat_f) / h_fg
    
    if x_eq <= 0:
        alpha = 0.0
    else:
        x_in = min(x_eq, 1.0)
        alpha = calculate_chexal_lellouche(x_in, G, P_sys)
        
    cases_data['Constant']['h'].append(h_z)
    cases_data['Constant']['x'].append(x_eq)
    cases_data['Constant']['alpha'].append(alpha)

# CASE 2: Sinusoidal Heat Flux
for z in z_axis:
    factor = (1 - np.cos(np.pi * z / L_rod)) / 2.0
    h_z_sin = h_inlet + (Q_channel / m_dot) * factor
    x_eq_sin = (h_z_sin - h_sat_f) / h_fg
    
    if x_eq_sin <= 0:
        alpha = 0.0
    else:
        x_in = min(x_eq_sin, 1.0)
        alpha = calculate_chexal_lellouche(x_in, G, P_sys)
    
    cases_data['Sinusoidal']['h'].append(h_z_sin)
    cases_data['Sinusoidal']['x'].append(x_eq_sin)
    cases_data['Sinusoidal']['alpha'].append(alpha)

# --- B. Pressure Drop Calculation Loop ---

pressure_results = {
    'Constant': {'grav': [], 'fric': [], 'acc': [], 'total': []},
    'Sinusoidal': {'grav': [], 'fric': [], 'acc': [], 'total': []}
}

for case_name in ['Constant', 'Sinusoidal']:
    
    # Init accumulators
    cum_grav = 0
    cum_fric = 0
    cum_acc = 0
    
    # Init Momentum Flux Tracker
    # Inlet assumption: x=0, alpha=0 -> Momentum flux = G^2/rho_f
    prev_mom_flux = (G**2) / rho_f 
    
    alphas = cases_data[case_name]['alpha']
    qualities = cases_data[case_name]['x']
    
    for i in range(len(z_axis)):
        alpha = alphas[i]
        x_eq = qualities[i]
        x_in = min(max(x_eq, 0.0), 1.0) # Clamp for property calc
        
        # 1. Gravity: dP = rho_bulk * g * dz
        rho_bulk = alpha * rho_g + (1 - alpha) * rho_f
        dP_g_step = rho_bulk * g_const * dz
        cum_grav += dP_g_step
        
        # 2. Friction: dP = Phi^2 * dP_lo
        # Liquid-only gradient
        Re_fo = (G * D_e) / mu_f
        f_fo = get_friction_factor_cheng_todreas_bundle(Re_fo)
        dP_lo_step = f_fo * (G**2) / (2 * rho_f * D_e) * dz
        
        if x_in <= 0:
            phi_sq = 1.0
        else:
            phi_sq = calculate_friedel_multiplier(x_in, G, D_e)
            
        dP_f_step = phi_sq * dP_lo_step
        cum_fric += dP_f_step
        
        # 3. Acceleration: dP = d(Momentum_Flux)
        # 1/rho_m = x^2/(rho_g*alpha) + (1-x)^2/(rho_f*(1-alpha))
        if x_in <= 0 or alpha <= 1e-6:
            curr_mom_term = 1.0 / rho_f
        else:
            curr_mom_term = (x_in**2)/(rho_g * alpha) + ((1-x_in)**2)/(rho_f * (1-alpha))
            
        curr_mom_flux = (G**2) * curr_mom_term
        
        # dP_acc is the change in momentum flux
        dP_a_step = curr_mom_flux - prev_mom_flux
        prev_mom_flux = curr_mom_flux
        cum_acc += dP_a_step
        
        # Store cumulative results
        pressure_results[case_name]['grav'].append(cum_grav)
        pressure_results[case_name]['fric'].append(cum_fric)
        pressure_results[case_name]['acc'].append(cum_acc)
        pressure_results[case_name]['total'].append(cum_grav + cum_fric + cum_acc)

# --- 5. Plotting (Void/Quality & Pressure) ---

# Plot 1: Void Fraction & Quality 
plt.figure(figsize=(10, 6))
plt.plot(z_axis, cases_data['Constant']['alpha'], 'b-', label='Void Fraction (Constant)')
plt.plot(z_axis, cases_data['Constant']['x'], 'b--', label='Quality (Constant)')
plt.plot(z_axis, cases_data['Sinusoidal']['alpha'], 'r-', label='Void Fraction (Sinusoidal)')
plt.plot(z_axis, cases_data['Sinusoidal']['x'], 'r--', label='Quality (Sinusoidal)')
plt.xlabel('Core Height (m)')
plt.ylabel('Fraction')
plt.title('Void Fraction & Quality Comparison')
plt.legend()
plt.grid(True)
_save_fig("void_fraction_quality")
plt.show()

# Plot 2: Pressure Losses Comparison


# Subplot 1: Constant Flux
plt.figure(figsize=(10, 6))
p_const = pressure_results['Constant']
plt.plot(z_axis, p_const['grav'], label='Gravity')
plt.plot(z_axis, p_const['fric'], label='Friction (Friedel)')
plt.plot(z_axis, p_const['acc'], label='Acceleration')
plt.plot(z_axis, p_const['total'], 'k--', linewidth=2, label='Total dP')
plt.title('Pressure Loss Components (Constant Flux)')
plt.xlabel('Height (m)')
plt.ylabel('Cumulative Pressure Drop (Pa)')
plt.legend()
plt.grid(True)
_save_fig("pressure_drop_constant")

# Subplot 2: Sinusoidal Flux
plt.figure(figsize=(10, 6))
p_sin = pressure_results['Sinusoidal']
plt.plot(z_axis, p_sin['grav'], label='Gravity')
plt.plot(z_axis, p_sin['fric'], label='Friction (Friedel)')
plt.plot(z_axis, p_sin['acc'], label='Acceleration')
plt.plot(z_axis, p_sin['total'], 'k--', linewidth=2, label='Total dP')
plt.title('Pressure Loss Components (Sinusoidal Flux)')
plt.xlabel('Height (m)')
plt.ylabel('Cumulative Pressure Drop (Pa)')
plt.legend()
plt.grid(True)
_save_fig("pressure_drop_sinusoidal")

plt.tight_layout()
plt.show()

# --- 6. Temperature & Saturation Height ---
T_sat = 285.8  # degrees C

def get_temperature_profile(h_profile, h_inlet, h_sat_f, T_inlet, T_sat):
    temps = []
    z_sat_found = None
    
    for i, h in enumerate(h_profile):
        if h < h_sat_f:
            fraction = (h - h_inlet) / (h_sat_f - h_inlet)
            T_val = T_inlet + fraction * (T_sat - T_inlet)
            temps.append(T_val)
        else:
            temps.append(T_sat)
            if z_sat_found is None:
                z_prev = z_axis[i-1]
                z_curr = z_axis[i]
                h_prev = h_profile[i-1]
                h_curr = h
                fraction_z = (h_sat_f - h_prev) / (h_curr - h_prev)
                z_sat_found = z_prev + fraction_z * (z_curr - z_prev)
    return temps, z_sat_found

temps_const, z_sat_const = get_temperature_profile(cases_data['Constant']['h'], h_inlet, h_sat_f, T_inlet, T_sat)
temps_sin, z_sat_sin = get_temperature_profile(cases_data['Sinusoidal']['h'], h_inlet, h_sat_f, T_inlet, T_sat)

# Plot Temps
plt.figure(figsize=(10, 6))
plt.plot(z_axis, temps_const, 'b-', label='Temp (Constant Flux)')
plt.plot(z_axis, temps_sin, 'r-', label='Temp (Sinusoidal Flux)')
plt.axhline(y=T_sat, color='k', linestyle=':', label='Saturation Temp')

if z_sat_const:
    plt.axvline(x=z_sat_const, color='b', linestyle='--', alpha=0.5)
    plt.text(z_sat_const, 280, f' Z_sat = {z_sat_const:.2f} m', color='b', rotation=90)
if z_sat_sin:
    plt.axvline(x=z_sat_sin, color='r', linestyle='--', alpha=0.5)
    plt.text(z_sat_sin, 280, f' Z_sat = {z_sat_sin:.2f} m', color='r', rotation=90)

plt.xlabel('Height (m)')
plt.ylabel('Temperature (°C)')
plt.title('Bulk Coolant Temperature vs Height')
plt.legend()
plt.grid(True)
_save_fig("temperature_profile")
plt.show()

print("-" * 30)
print(f"Saturation Height (Constant):   {z_sat_const:.4f} m")
print(f"Saturation Height (Sinusoidal): {z_sat_sin:.4f} m")
print("-" * 30)
print(f"Total Pressure Drop (Constant):   {pressure_results['Constant']['total'][-1]/1e5:.3f} bar")
print(f"Total Pressure Drop (Sinusoidal): {pressure_results['Sinusoidal']['total'][-1]/1e5:.3f} bar")
print("-" * 30)

# --- 7. Dryout Power Search ---
print("7: DRYOUT SEARCH")
# Constants
current_power = 4500e6  # Start at nominal power
power_step = 100e6      # 100 MWt increment [cite: 21]

# Lists for Plotting
pin_powers = []       # X-axis: Average Pin Power (kW)
exit_voids = []       # Y-axis 1: Exit Void Fraction
exit_qualities = []   # Y-axis 2: Exit Quality

dryout_found = False

while not dryout_found:
    # 1. Calculate Power Parameters
    Q_channel = current_power / N_rods
    
    # 2. Calculate Exit Enthalpy (Sinusoidal Profile is most limiting/realistic)
    # At exit (z = L), the cosine integral factor is exactly 1.0:
    # Factor = (1 - cos(pi * L/L))/2 = (1 - (-1))/2 = 1.0
    # So h_exit = h_in + Q_channel / m_dot
    h_exit = h_inlet + (Q_channel / m_dot)
    
    # 3. Calculate Exit Conditions
    x_exit = (h_exit - h_sat_f) / h_fg
    
    # Clamp for void calculation (so code doesn't crash if x > 1)
    x_in = min(max(x_exit, 0.0), 1.0)
    
    if x_in <= 0:
        alpha_exit = 0.0
    else:
        alpha_exit = calculate_chexal_lellouche(x_in, G, P_sys)
        
    # 4. Store Data
    # Convert Pin Power to kW for cleaner plotting
    pin_powers.append(Q_channel / 1000.0) 
    exit_qualities.append(x_exit)
    exit_voids.append(alpha_exit)
    
    # 5. Check for Dryout Condition
    if x_exit >= 1.0:
        dryout_found = True
        print(f"Core Thermal Power: {current_power/1e6:.0f} MWt")
        print(f"Average Pin Power:  {Q_channel/1000.0:.2f} kW")
        print(f"Exit Quality:       {x_exit:.4f}")
    else:
        # Increment and continue
        current_power += power_step
        
        # Safety break to prevent infinite loop if physics fail
        if current_power > 100000e6: 
            print("Loop stopped: Power exceeded 100,000 MWt without dryout.")
            break

# --- Plotting Task 7 ---
fig, ax1 = plt.subplots(figsize=(10, 6))

color = 'tab:red'
ax1.set_xlabel('Average Pin Power (kW)')
ax1.set_ylabel('Exit Quality', color=color)
ax1.plot(pin_powers, exit_qualities, color=color, linewidth=2, label='Exit Quality')
ax1.tick_params(axis='y', labelcolor=color)
ax1.axhline(y=1.0, color='k', linestyle=':', label='Dryout Limit (x=1)')

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
color = 'tab:blue'
ax2.set_ylabel('Exit Void Fraction', color=color)  # we already handled the x-label with ax1
ax2.plot(pin_powers, exit_voids, color=color, linewidth=2, linestyle='--', label='Exit Void Fraction')
ax2.tick_params(axis='y', labelcolor=color)

plt.title('Exit Conditions vs. Rod Power (Dryout Search)')
fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.grid(True)
_save_fig("dryout_search")
plt.show()