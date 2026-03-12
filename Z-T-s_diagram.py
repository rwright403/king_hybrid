import numpy as np
import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP

fluid = "N2O"

# -------------------------
# Saturation dome
# -------------------------
T_trip = CP.PropsSI("Ttriple", fluid)
T_crit = CP.PropsSI("Tcrit", fluid)

T_sat = np.linspace(T_trip + 0.5, T_crit - 0.5, 400)

s_f, s_g, T_dome = [], [], []
for T in T_sat:
    try:
        s_l = CP.PropsSI("S", "T", T, "Q", 0, fluid)
        s_v = CP.PropsSI("S", "T", T, "Q", 1, fluid)
        s_f.append(s_l)
        s_g.append(s_v)
        T_dome.append(T)
    except:
        pass

# -------------------------
# Isobars (constant P)
# -------------------------
isobars_MPa = [0.3, 1.3, 3.5, 6.0, 8.3]   # MPa
T_grid = np.linspace(200, 400, 400)       # K

# -------------------------
# Isenthalps (constant h)
# NOTE: choose values that actually intersect your plotted region
# -------------------------
isenthalps_kJkg = [50, 150, 250, 350, 450]   # kJ/kg

# Pressure sweep for isenthalps
P_min = 0.1e6    # Pa
P_max = 9.5e6    # Pa
P_grid = np.linspace(P_min, P_max, 600)

# -------------------------
# Plot
# -------------------------
plt.figure(figsize=(10, 7))

# Saturation dome
plt.plot(s_f, T_dome, "k", lw=2)
plt.plot(s_g, T_dome, "k", lw=2)

# Isobars
for P_MPa in isobars_MPa:
    P = P_MPa * 1e6
    s_vals, T_vals = [], []
    for T in T_grid:
        try:
            s = CP.PropsSI("S", "T", T, "P", P, fluid)
            s_vals.append(s)
            T_vals.append(T)
        except:
            pass
    if len(s_vals) > 10:
        plt.plot(s_vals, T_vals, color="red", lw=1.5)
        mid = len(s_vals)//2
        plt.text(s_vals[mid], T_vals[mid], f"P = {P_MPa:.1f} MPa", color="red")

# Isenthalps (TRUE constant h: sweep P, compute T(P), s(P) at fixed h)
for h_kJkg in isenthalps_kJkg:
    h = h_kJkg * 1e3  # J/kg
    s_vals, T_vals = [], []
    for P in P_grid:
        try:
            T = CP.PropsSI("T", "P", P, "H", h, fluid)
            s = CP.PropsSI("S", "P", P, "H", h, fluid)

            # Keep only points in your plotted window (helps avoid weird far-out values)
            if 200 <= T <= 400 and 0 <= s <= 3000:
                s_vals.append(s)
                T_vals.append(T)
        except:
            pass

    if len(s_vals) > 10:
        # sort by entropy so the line draws nicely
        idx = np.argsort(s_vals)
        s_vals = np.array(s_vals)[idx]
        T_vals = np.array(T_vals)[idx]

        plt.plot(s_vals, T_vals, color="magenta", lw=2)

        # label near the top-ish point if possible
        k = np.argmax(T_vals)
        plt.text(s_vals[k], T_vals[k], f"h = {h_kJkg:.1f} kJ/kg", color="magenta")

# Formatting
plt.xlabel("Entropy (J/kg·K)")
plt.ylabel("Temperature (K)")
plt.title("T–s Diagram for Nitrous Oxide (N₂O) with Isobars and Isenthalps")
plt.grid(True)
plt.tight_layout()
plt.show()
