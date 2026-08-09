# Overall-Heat-Transfer-Coefficeint-Calcilator
Overall heat transfer coefficient (U) for tubes embedded in a sand packed-bed TES. Combines Gnielinski/Dittus–Boelter (single-phase), Gungor–Winterton + Cooper (flow boiling), radial wall conduction (Incropera), and stagnant-bed sand-side resistance (Nu = 2).

"""
Overall heat transfer coefficient U for embedded tubes in a sand packed-bed TES.

Resistance network (per unit tube length):
    1/(U*A) = R_i + R_wall + R_sand
            = 1/(h_i * π * Di) + ln(Do/Di)/(2π k_cu) + 1/(h_sand * π * Do)

Inner-side h_i is calculated with:
  - Single-phase (water or steam): Gnielinski / Dittus-Boelter
  - Flow boiling: Gungor-Winterton (1986) + Cooper pool-boiling term

Requires:  pip install CoolProp numpy pandas
"""

import numpy as np
import pandas as pd
from CoolProp.CoolProp import PropsSI

# ======================================================================
# Geometry / materials / operating conditions
# ======================================================================
Di, Do   = 0.070, 0.078          # m  (70 mm ID, 4 mm wall → Do = 78 mm)
K_CU     = 400.0                 # W/m·K  copper
K_EFF    = 1.5                   # W/m·K  effective conductivity of dry sand bed
MDOT, P  = 0.1, 1.0e6            # kg/s, Pa  (10 bar)
FLUID    = "Water"
M_WATER  = 18.015                # kg/kmol (for Cooper correlation)


# ======================================================================
# Fluid property helper
# ======================================================================
def props(T=None, p=None, Q=None):
    """Return dict of properties. Use (T,p) for single-phase or (p,Q) for saturation."""
    out = {}
    if Q is None:
        out["rho"] = PropsSI("D", "T", T, "P", p, FLUID)
        out["mu"]  = PropsSI("V", "T", T, "P", p, FLUID)
        out["k"]   = PropsSI("L", "T", T, "P", p, FLUID)
        out["cp"]  = PropsSI("C", "T", T, "P", p, FLUID)
    else:
        out["rho"] = PropsSI("D", "P", p, "Q", Q, FLUID)
        out["mu"]  = PropsSI("V", "P", p, "Q", Q, FLUID)
        out["k"]   = PropsSI("L", "P", p, "Q", Q, FLUID)
        out["cp"]  = PropsSI("C", "P", p, "Q", Q, FLUID)
    out["Pr"] = out["mu"] * out["cp"] / out["k"]
    return out


# ======================================================================
# Inner heat-transfer coefficients
# ======================================================================
def Re_D(mdot, D, mu):
    return 4.0 * mdot / (np.pi * D * mu)


def dittus_boelter(Re, Pr, heating=True):
    n = 0.4 if heating else 0.3
    return 0.023 * Re**0.8 * Pr**n


def gnielinski(Re, Pr):
    f = (0.79 * np.log(Re) - 1.64)**-2
    return (f/8) * (Re - 1000) * Pr / (1 + 12.7*(f/8)**0.5 * (Pr**(2/3) - 1))


def h_single_phase(mdot, D, T, p, correlation="gnielinski"):
    """h [W/m²K] for single-phase liquid water or superheated steam."""
    pr = props(T=T, p=p)
    Re = Re_D(mdot, D, pr["mu"])
    Nu = gnielinski(Re, pr["Pr"]) if correlation == "gnielinski" else dittus_boelter(Re, pr["Pr"])
    return Nu * pr["k"] / D, Re, pr["Pr"]


def h_gungor_winterton(G, D, p, x, qpp):
    """
    Two-phase flow boiling h [W/m²K] – Gungor-Winterton (1986) + Cooper (1984).
    G   : mass flux [kg/m²s]
    D   : inner diameter [m]
    p   : pressure [Pa]
    x   : vapor quality [-]
    qpp : wall heat flux [W/m²]
    """
    liq  = props(p=p, Q=0)
    vap  = props(p=p, Q=1)
    h_fg = PropsSI("H", "P", p, "Q", 1, FLUID) - PropsSI("H", "P", p, "Q", 0, FLUID)
    pc   = PropsSI("pcrit", FLUID)

    Re_l = G * (1 - x) * D / liq["mu"]
    h_L  = 0.023 * Re_l**0.8 * liq["Pr"]**0.4 * liq["k"] / D

    Bo   = qpp / (G * h_fg)
    Xtt  = ((1-x)/x)**0.9 * (vap["rho"]/liq["rho"])**0.5 * (liq["mu"]/vap["mu"])**0.1
    E    = 1 + 24000*Bo**1.16 + 1.37*(1/Xtt)**0.86
    S    = 1.0 / (1.0 + 1.15e-6 * E**2 * Re_l**1.17)

    pr_red = p / pc
    h_nb = 55 * pr_red**0.12 * (-np.log10(pr_red))**-0.55 * M_WATER**-0.5 * qpp**0.67

    h_tp = E * h_L + S * h_nb
    return h_tp, dict(h_L=h_L, E=E, S=S, h_nb=h_nb, Bo=Bo, Xtt=Xtt, Re_l=Re_l)


# ======================================================================
# Overall U (resistance network)
# ======================================================================
def R_wall(Do, Di, k_wall=K_CU):
    """Radial wall conduction resistance per unit length [m·K/W]."""
    return np.log(Do / Di) / (2 * np.pi * k_wall)


def h_sand(Do, k_eff=K_EFF):
    """Sand-side coefficient – stagnant conduction limit (Nu = 2)."""
    return 2.0 * k_eff / Do


def overall_U(h_i, Do=Do, Di=Di, k_wall=K_CU, k_eff=K_EFF, basis="outer"):
    """
    Return U [W/m²K] on outer- or inner-area basis together with resistance shares.
    """
    hs = h_sand(Do, k_eff)
    Ri = 1.0 / (h_i * np.pi * Di)
    Rw = R_wall(Do, Di, k_wall)
    Rs = 1.0 / (hs * np.pi * Do)
    Rt = Ri + Rw + Rs
    A  = np.pi * (Do if basis == "outer" else Di)
    U  = 1.0 / (Rt * A)
    return U, dict(Ri=Ri, Rw=Rw, Rs=Rs, h_sand=hs,
                   share_inner=Ri/Rt, share_wall=Rw/Rt, share_sand=Rs/Rt)


# ======================================================================
# Main evaluation
# ======================================================================
if __name__ == "__main__":
    A = np.pi * Di**2 / 4
    G = MDOT / A

    regimes = {}

    # 1. Sensible water heating (60 °C feedwater)
    h, _, _ = h_single_phase(MDOT, Di, T=333.15, p=P)
    regimes["1. Sensible water (60°C)"] = h

    # 2. Flow boiling at x = 0.5, q'' = 3.66 kW/m²
    h, _ = h_gungor_winterton(G, Di, P, x=0.5, qpp=3.66e3)
    regimes["2. Flow boiling (x=0.5)"] = h

    # 3. Superheated steam (500 °C, 10 bar)
    h, _, _ = h_single_phase(MDOT, Di, T=773.15, p=P)
    regimes["3. Superheated steam (500°C)"] = h

    # ---- Table of U for each regime ----
    rows = []
    for name, hi in regimes.items():
        U, info = overall_U(hi)
        rows.append(dict(
            regime=name,
            h_i=round(hi),
            h_sand=round(info["h_sand"], 1),
            U_outer=round(U, 1),
            sand_share_pct=round(100 * info["share_sand"], 1)
        ))
    df = pd.DataFrame(rows)
    print(df.to_string(index=False))

    # ---- Sensitivity: U vs effective sand conductivity (boiling regime) ----
    print("\nSensitivity of U (boiling regime) to sand k_eff:")
    for ke in [0.25, 0.4, 0.6, 1.0, 2.0, 4.0]:
        U, info = overall_U(regimes["2. Flow boiling (x=0.5)"], k_eff=ke)
        print(f"  k_eff = {ke:4.2f} W/m·K  →  h_sand = {info['h_sand']:6.1f},  U = {U:6.1f} W/m²K")

    df.to_csv("overall_U_results.csv", index=False)
    print("\nSaved → overall_U_results.csv")
